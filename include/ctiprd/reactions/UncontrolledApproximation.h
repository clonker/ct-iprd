/**
 * @file UncontrolledApproximation.h
 * @brief 
 * @author clonker
 * @date 3/27/22
 */
#pragma once

#include <ctiprd/any.h>
#include <ctiprd/NeighborList.h>
#include "basic.h"

namespace ctiprd::reactions {

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

        auto worker = [&system](auto id, typename Particles::Position &pos, const typename Particles::ParticleType &type,
                         typename Particles::Force &force) {
            std::apply([&pos, &type, &force](auto &&... reaction) {
                (
                        (force += reaction.force(pos, type)),
                        ...);
            }, system.reactionsO1);
        };

        particles->forEachParticle(worker, pool);
        pool->waitForTasks();
    }

    std::unique_ptr<NeighborList> neighborList_;
};

}
