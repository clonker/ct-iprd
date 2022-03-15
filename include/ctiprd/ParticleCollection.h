//
// Created by mho on 5/5/21.
//

#pragma once

#include <cstdint>
#include <vector>
#include <map>
#include <future>
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
    using Force = Vec<dtype, DIM>;
    using Velocity = Vec<dtype, DIM>;

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

        auto n = nParticles();
        auto *pos = positions_.data();
        auto *forces = forces_.data();
        auto *velocities = velocities_.data();

        Vec<dtype, DIM> defaultPosition {};
        Vec<dtype, DIM> defaultForce {};
        Vec<dtype, DIM> defaultVelocity {};

        auto grainSize = n / granularity;

        auto loop = [operation = std::forward<F>(op), pos, forces, velocities](std::size_t beginIx, std::size_t endIx) {
            for (std::size_t i = beginIx; i < endIx; ++i) {
                if constexpr(containsForces() && containsVelocities()) {
                    operation(i, *(pos + i), *(forces + i), *(velocities + i));
                } else if constexpr(containsForces() && !containsVelocities()) {
                    operation(i, *(pos + i), *(forces + i));
                } else if constexpr(!containsForces() && containsVelocities()) {
                    operation(i, *(pos + i), *(velocities + i));
                } else {
                    operation(i, *(pos + i));
                }
            }
        };
        std::size_t beginIx = 0;
        for (std::size_t i = 0; i < granularity - 1; ++i) {
            auto endIx = beginIx + grainSize;
            if (beginIx != endIx) {
                futures.push_back(pool->push(loop, beginIx, endIx));
            }
            beginIx = endIx;
        }
        if (beginIx != n) {
            futures.push_back(pool->push(loop, beginIx, n));
        }

        return std::move(futures);
    }

    [[nodiscard]] std::string_view particleType() const {
        return particleType_;
    }

    const std::vector<Position> &positions () const {
        return positions_;
    }

    const std::vector<Force> &forces() const {
        return forces_;
    }

private:
    std::vector<Position> positions_;
    std::vector<Force> forces_;
    std::vector<Velocity> velocities_;
    std::string_view particleType_;
};
}
