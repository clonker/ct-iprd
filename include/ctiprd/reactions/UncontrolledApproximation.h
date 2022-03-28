/**
 * @file UncontrolledApproximation.h
 * @brief 
 * @author clonker
 * @date 3/27/22
 */
#pragma once

#include <ctiprd/any.h>
#include <ctiprd/NeighborList.h>
#include <variant>
#include "basic.h"

namespace ctiprd::reactions {

template<typename dtype>
struct ReactionEvent {
    std::size_t id1, id2;
    std::variant<const doi::ReactionE1<dtype>*, const doi::ReactionE2<dtype>*> reaction;
};

template<typename System>
struct UncontrolledApproximation {
    static constexpr std::size_t DIM = System::DIM;
    using dtype = typename System::dtype;
    using NeighborList = ctiprd::nl::NeighborList<DIM, System::periodic, dtype>;

    using ReactionsO1 = typename System::ReactionsO1;
    using ReactionsO2 = typename System::ReactionsO2;

    static constexpr int nReactionsO1 = std::tuple_size_v<ReactionsO1>;
    static constexpr int nReactionsO2 = std::tuple_size_v<ReactionsO2>;

    template<typename Particles, typename Pool>
    void reactions(const System &system, dtype tau, std::shared_ptr<Particles> particles,
                   std::shared_ptr<Pool> pool) {
        if constexpr(nReactionsO2 > 0) {
            if (!neighborList_) {
                auto cutoff = reactions::reactionRadius<dtype>(system.reactionsO2);
                neighborList_ = std::make_unique<NeighborList>(System::boxSize, cutoff);
            }
            neighborList_->update(particles, pool);
        }

        auto worker = [&system, tau](auto id, typename Particles::Position &pos, const typename Particles::ParticleType &type,
                         typename Particles::Force &force) {
            std::vector<ReactionEvent<dtype>> events;

            std::apply([&pos, &type, &events, &tau, &id](auto &&... reaction) {
                (
                        ([&events, &pos, &type, &tau, &id](const auto & r) {
                            if(r.shouldPerform(tau, pos, type)) {
                                events.push_back(ReactionEvent<dtype>{
                                    .id1 = id,
                                    .id2 = id,
                                    .reaction = &r
                                });
                            }
                        }(reaction)),
                        ...);
            }, system.reactionsO1);

            return events;
        };

        std::vector<ReactionEvent<dtype>> events;
        auto futures = particles->forEachParticle(worker, pool);
        for(auto &future : futures) {
            auto evts = future.get();
            for (const auto &x : evts) {
                events.reserve(events.size() + evts.size());
                events.insert(end(events), begin(x), end(x));
            }
        }
    }

    std::unique_ptr<NeighborList> neighborList_;
};

}
