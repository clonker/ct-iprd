/**
 *
 *
 * @file ForceField.h
 * @brief 
 * @author clonker
 * @date 3/27/22
 */
#pragma once

#include <tuple>

namespace ctiprd::potentials {

template<typename System>
struct ForceField {
    static constexpr std::size_t DIM = System::DIM;
    using dtype = typename System::dtype;
    using ExternalPotentials = typename System::ExternalPotentials;
    using PairPotentials = typename System::PairPotentials;

    static constexpr int nExternalPotentials = std::tuple_size_v<ExternalPotentials>;
    static constexpr int nPairPotentials = std::tuple_size_v<PairPotentials>;
    using NeighborList = nl::NeighborList<DIM, System::periodic, dtype>;

    void setExternalPotentials(const ExternalPotentials &potentials) {
        externalPotentials_ = potentials;
    }

    void setPairPotentials(const PairPotentials &potentials) {
        pairPotentials_ = potentials;
        neighborList_.reset();
    }

    template<typename Particles, typename Pool>
    void forces(std::shared_ptr<Particles> particles, std::shared_ptr<Pool> pool, bool wait = true) {
        if constexpr(nPairPotentials > 0) {
            if (!neighborList_) {
                auto cutoff = potentials::cutoff<dtype>(pairPotentials_);
                neighborList_ = std::make_unique<NeighborList>(System::boxSize, cutoff);
            }
            neighborList_->update(particles, pool);
        }

        auto worker = [
                &pot = externalPotentials_,
                &potPair = pairPotentials_,
                nl = neighborList_.get(),
                &data = *particles
        ]
                (auto id, typename Particles::Position &pos, const typename Particles::ParticleType &type,
                 typename Particles::Force &force) {

            std::fill(begin(force.data), end(force.data), static_cast<dtype>(0));
            std::apply([&pos, &type, &force](auto &&... args) {
                ((force += args.force(pos, type)), ...);
            }, pot);

            if constexpr(nPairPotentials > 0) {
                nl->forEachNeighbor(id, data, [&pos, &type, &force, &potPair](auto neighborId, const auto &neighborPos,
                                                                              const auto &neighborType,
                                                                              const auto &neighborForce) {
                    std::apply([&pos, &type, &force, &neighborPos, &neighborType](auto &&... args) {
                        ((force += args.force(pos, type, neighborPos, neighborType)), ...);
                    }, potPair);
                });
            }
        };
        particles->forEachParticle(worker, pool);
        if (wait) {
            pool->waitForTasks();
        }
    }

    ExternalPotentials externalPotentials_;
    PairPotentials pairPotentials_;
    std::unique_ptr<NeighborList> neighborList_;

};

}
