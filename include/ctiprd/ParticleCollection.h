//
// Created by mho on 5/5/21.
//

#pragma once

#include <cstdint>
#include <vector>
#include <map>
#include <future>
#include <optional>
#include <bitset>

#include <ctiprd/vec.h>
#include <ctiprd/config.h>
#include "ctiprd/systems/util.h"

namespace ctiprd {

namespace particle_collection {
enum {
    usePositions = 0b001,
    useForces = 0b010,
    useVelocities = 0b100,
};
}

template<typename System, int flags = particle_collection::usePositions | particle_collection::useForces>
class ParticleCollection {
public:
    using dtype = typename System::dtype;
    static constexpr std::size_t DIM = System::DIM;

    using Position = Vec<dtype, DIM>;
    using MaybePosition = std::optional<Position>;
    using Force = Vec<dtype, DIM>;
    using Velocity = Vec<dtype, DIM>;
    using ParticleType = std::size_t;

    template<typename T>
    using ContainerType = std::deque<T>;

    using size_type = typename ContainerType<MaybePosition>::size_type;

    static constexpr int dim = DIM;

    static constexpr bool containsPositions() { return flags & particle_collection::usePositions; }
    static constexpr bool containsForces() { return flags & particle_collection::useForces; }
    static constexpr bool containsVelocities() { return flags & particle_collection::useVelocities; }

    [[nodiscard]] std::size_t nParticles() const {
        return positions_.size();
    }

    void addParticle(const Position &position, const char* type) {
        positions_.push_back(position);
        if constexpr(containsForces()) {
            forces_.template emplace_back();
        }
        if constexpr(containsVelocities()) {
            velocities_.template emplace_back();
        }
        particleTypes_.push_back(systems::particleTypeId<System::types>(type));
    }

    [[nodiscard]] ParticleType typeOf(std::size_t ix) const {
        return particleTypes_[ix];
    }

    template<typename F, typename Pool>
    std::vector<std::future<void>> forEachParticle(F &&op, config::PoolPtr<Pool> pool) {
        std::vector<std::future<void>> futures;
        auto granularity = config::threadGranularity(pool);
        futures.reserve(granularity);

        auto loop = [operation = std::forward<F>(op)](
                auto startIndex,
                const auto &beginPositions, const auto &endPositions,
                auto itTypes, auto itForces, auto itVelocities
        ) {
            for(auto itPos = beginPositions; itPos != endPositions; ++itPos, ++startIndex) {
                if(*itPos) {
                    if constexpr(containsForces() && containsVelocities()) {
                        operation(startIndex, **itPos, *itTypes, *itForces, *itVelocities);
                    } else if constexpr(containsForces() && !containsVelocities()) {
                        operation(startIndex, **itPos, *itTypes, *itForces);
                    } else if constexpr(!containsForces() && containsVelocities()) {
                        operation(startIndex, **itPos, *itTypes, *itVelocities);
                    } else {
                        operation(startIndex, **itPos, *itTypes);
                    }
                }

                if constexpr(containsForces()) {
                    ++itForces;
                }
                if constexpr(containsVelocities()) {
                    ++itVelocities;
                }
            }
        };

        auto n = nParticles();
        auto grainSize = n / granularity;
        auto itPos = positions_.begin();
        auto itTypes = particleTypes_.begin();
        auto itForces = forces_.begin();
        auto itVelocities = velocities_.begin();

        std::size_t beginIx = 0;
        for (std::size_t i = 0; i < granularity - 1; ++i) {
            if (itPos != positions_.end()) {
                futures.push_back(pool->push(loop, beginIx, itPos, itPos + grainSize, itTypes, itForces, itVelocities));
            }
            beginIx += grainSize;
            itPos += grainSize;
            itTypes += grainSize;
            if constexpr(containsForces()) {
                itForces += grainSize;
            }
            if constexpr(containsVelocities()) {
                itVelocities += grainSize;
            }
        }
        if (itPos != positions_.end()) {
            futures.push_back(pool->push(loop, beginIx, itPos, positions_.end(), itTypes, itForces, itVelocities));
        }

        return std::move(futures);
    }

    const ContainerType<MaybePosition> &positions () const {
        return positions_;
    }

    const ContainerType<Force> &forces() const {
        return forces_;
    }

    const Position &position(size_type index) const {
        return *positions_[index];
    }

    const Force &force(size_type index) const {
        return forces_[index];
    }

private:
    ContainerType<MaybePosition> positions_;
    ContainerType<Force> forces_;
    ContainerType<Velocity> velocities_;
    ContainerType<ParticleType> particleTypes_;
};
}
