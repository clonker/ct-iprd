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

namespace ctiprd {

using ParticleType = std::uint32_t;

namespace particle_collection {
enum {
    usePositions = 0b001,
    useForces = 0b010,
    useVelocities = 0b100,
};
}

template<int DIM, typename dtype, int flags = particle_collection::usePositions | particle_collection::useForces>
class ParticleCollection {
public:
    using Position = Vec<dtype, DIM>;
    using MaybePosition = std::optional<Position>;
    using Force = Vec<dtype, DIM>;
    using Velocity = Vec<dtype, DIM>;

    template<typename T>
    using ContainerType = std::deque<T>;

    static constexpr int dim = DIM;

    static constexpr bool containsPositions() { return flags & particle_collection::usePositions; }
    static constexpr bool containsForces() { return flags & particle_collection::useForces; }
    static constexpr bool containsVelocities() { return flags & particle_collection::useVelocities; }

    explicit ParticleCollection(std::string_view particleType) : particleType_(particleType) {}

    [[nodiscard]] std::size_t nParticles() const {
        return positions_.size();
    }

    void addParticle(const Position &position) {
        positions_.push_back(position);
        forces_.template emplace_back();
        velocities_.template emplace_back();
    }

    template<typename F, typename Pool>
    std::vector<std::future<void>> forEachParticle(F &&op, config::PoolPtr<Pool> pool) {
        std::vector<std::future<void>> futures;
        auto granularity = config::threadGranularity(pool);
        futures.reserve(granularity);


        Vec<dtype, DIM> defaultPosition {};
        Vec<dtype, DIM> defaultForce {};
        Vec<dtype, DIM> defaultVelocity {};


        auto loop = [operation = std::forward<F>(op)](
                auto startIndex,
                auto beginPositions, auto endPositions,
                auto itForces, auto itVelocities
        ) {
            for(auto itPos = beginPositions; itPos != endPositions; ++itPos, ++startIndex) {
                if constexpr(containsForces() && containsVelocities()) {
                    operation(startIndex, *itPos, *itForces, *itVelocities);
                } else if constexpr(containsForces() && !containsVelocities()) {
                    operation(startIndex, *itPos, *itForces);
                } else if constexpr(!containsForces() && containsVelocities()) {
                    operation(startIndex, *itPos, *itVelocities);
                } else {
                    operation(startIndex, *itPos);
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
        auto itForces = forces_.begin();
        auto itVelocities = velocities_.begin();

        std::size_t beginIx = 0;
        for (std::size_t i = 0; i < granularity - 1; ++i) {
            if (itPos != positions_.end()) {
                futures.push_back(pool->push(loop, beginIx, itPos, itPos + grainSize, itForces, itVelocities));
            }
            beginIx += grainSize;
            itPos += grainSize;
            if constexpr(containsForces()) {
                itForces += grainSize;
            }
            if constexpr(containsVelocities()) {
                itVelocities += grainSize;
            }
        }
        if (itPos != positions_.end()) {
            futures.push_back(pool->push(loop, beginIx, itPos, positions_.end(), itForces, itVelocities));
        }

        return std::move(futures);
    }

    [[nodiscard]] std::string_view particleType() const {
        return particleType_;
    }

    const ContainerType<MaybePosition> &positions () const {
        return positions_;
    }

    const ContainerType<Force> &forces() const {
        return forces_;
    }

private:
    ContainerType<MaybePosition> positions_;
    ContainerType<Force> forces_;
    ContainerType<Velocity> velocities_;
    std::string_view particleType_;
};
}
