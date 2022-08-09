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

#include <ctiprd/potentials/util.h>
#include <ctiprd/potentials/external.h>
#include <ctiprd/potentials/interaction.h>

#include <ctiprd/cpu/NeighborList.h>

#include "reactions.h"
#include "forces.h"

namespace ctiprd::cpu::potentials {

template<typename ParticleCollection, systems::system System>
struct ForceField {
    static constexpr std::size_t DIM = System::DIM;
    using dtype = typename System::dtype;
    using ExternalPotentials = typename System::ExternalPotentials;
    using PairPotentials = typename System::PairPotentials;

    static constexpr int nExternalPotentials = std::tuple_size_v<ExternalPotentials>;
    static constexpr int nPairPotentials = std::tuple_size_v<PairPotentials>;
    using NeighborList = nl::NeighborList<DIM, System::periodic, dtype>;

    template<typename Particles, typename Pool>
    void forces(std::shared_ptr<Particles> particles, std::shared_ptr<Pool> pool, bool wait = true) {
        if constexpr(nPairPotentials > 0) {
            neighborList_->update(particles.get(), pool);
        }

        if constexpr(nExternalPotentials > 0 || nPairPotentials > 0) {
            const auto worker = [
                    &pot = potentialsO1,
                    &potPair = potentialsO2,
                    nl = neighborList_.get(),
                    &data = *particles
            ]
                    (const auto &particleId, typename Particles::Position &pos, const typename Particles::ParticleType &type,
                     typename Particles::Force &force) {

                std::fill(begin(force.data), end(force.data), static_cast<dtype>(0));
                for(const auto &potentialO1 : pot[type]) {
                    force += potentialO1->force(pos);
                }

                if constexpr(nPairPotentials > 0) {
                    if(nl->isAllowedType(type)) {
                        nl->forEachNeighbor(particleId, data, [&pos, &type, &force, &potPair](auto neighborId, const auto &neighborPos,
                                                                                              const auto &neighborType,
                                                                                              const auto &neighborForce) {
                            for(const auto &potentialO2 : potPair[std::tie(type, neighborType)]) {
                                force += potentialO2->force(pos, neighborPos);
                            }
                        });
                    }
                }
            };
            auto futures = particles->forEachParticle(worker, pool);
            if (wait) {
                for(auto &future : futures)  {
                    future.wait();
                }
            }
        }
    }

    explicit ForceField(const System &system) {
        std::tie(potentialsO1, backingO1) = forces::generateMapO1<ParticleCollection>(system);
        std::tie(potentialsO2, backingO2) = forces::generateMapO2<ParticleCollection>(system);

        if constexpr(nPairPotentials > 0) {
            auto cutoff = ctiprd::potentials::cutoff<dtype>(system.pairPotentials);
            neighborList_ = std::make_unique<NeighborList>(System::boxSize, cutoff);

            std::unordered_set<std::size_t> neighborListTypes{};
            for (const auto &reactionElement: potentialsO2) {
                if (!reactionElement.second.empty()) {
                    neighborListTypes.emplace(std::get<0>(reactionElement.first));
                    neighborListTypes.emplace(std::get<1>(reactionElement.first));
                }
            }
            neighborList_->setTypes(neighborListTypes);
        }
    }

    forces::FFO1Map<ParticleCollection> potentialsO1;
    forces::FFO1Backing<ParticleCollection> backingO1;
    forces::FFO2Map<ParticleCollection> potentialsO2;
    forces::FFO2Backing<ParticleCollection> backingO2;

    std::unique_ptr<NeighborList> neighborList_;

};

}
