/**
 * @file UncontrolledApproximation.h
 * @brief 
 * @author clonker
 * @date 3/27/22
 */
#pragma once

#include <variant>
#include <execution>
#include <algorithm>

#include <ctiprd/any.h>
#include <ctiprd/NeighborList.h>
#include "basic.h"

namespace ctiprd::reactions {

template<typename dtype>
struct ReactionEvent {
    std::size_t id1, id2;
    std::variant<const doi::ReactionE1<dtype>*, const doi::ReactionE2<dtype>*> reaction;

    bool conflict(const ReactionEvent &other) const {
        return id1 == other.id1 || id1 == other.id2 || id2 == other.id1 || id2 == other.id2;
    }
};

template<typename System, typename Generator>
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

        std::mutex mutex;
        std::vector<ReactionEvent<dtype>> events;
        auto worker = [&system, tau, nl=neighborList_.get(), data=particles.get(), &events, &mutex](
                auto id, typename Particles::Position &pos, const typename Particles::ParticleType &type,
                typename Particles::Force &force
        ) {
            std::vector<ReactionEvent<dtype>> localEvents;

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

            if constexpr(nReactionsO2 > 0) {
                nl->template forEachNeighbor(id, *data, [&](auto neighborId, const auto &neighborPos,
                                                           const auto &neighborType,
                                                           const auto &neighborForce) {
                    std::apply([&](auto &&... reaction) {
                        (
                                ([&](const auto & r) {
                                    if(r.shouldPerform(tau, pos, type, neighborPos, neighborType)) {
                                        events.push_back(ReactionEvent<dtype>{
                                                .id1 = id,
                                                .id2 = neighborId,
                                                .reaction = &r
                                        });
                                    }
                                }(reaction)),
                                ...);
                    }, system.reactionsO2);
                });
            }

            {
                std::scoped_lock lock {mutex};
                events.reserve(events.size() + localEvents.size());
                events.insert(end(events), begin(localEvents), end(localEvents));
            }
        };

        particles->forEachParticle(worker, pool);
        pool->waitForTasks();
        std::shuffle(begin(events), end(events), rnd::staticThreadLocalGenerator<Generator>());
        ParticleCollectionUpdater updater {*particles};
        for (auto it = begin(events); it != end(events); ++it) {

        }
        particles->update(begin(updater.toAdd), end(updater.toAdd), begin(updater.toRemove), end(updater.toRemove));
    }

    std::unique_ptr<NeighborList> neighborList_;
};

}
