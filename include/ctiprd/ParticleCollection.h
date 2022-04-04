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

    [[nodiscard]] size_type size() const {
        return positions_.size();
    }

    void setType(size_type index, const ParticleType &type) {
        particleTypes_[index] = type;
    }

    void setPosition(size_type index, const Position &position) {
        positions_[index] = position;
    }

    template<typename T>
    void addParticle(const Position &position, T &&type) requires std::convertible_to<T, std::string_view> {
        addParticle(position, systems::particleTypeId<System::types>(type));
    }

    void addParticle(const Position &position, std::size_t type) {
        if (blanks.empty()) {
            positions_.push_back(position);
            if constexpr(containsForces()) {
                forces_.emplace_back();
            }
            if constexpr(containsVelocities()) {
                velocities_.emplace_back();
            }
            particleTypes_.push_back(type);
        } else {
            auto ix = blanks.back();
            positions_[ix] = position;
            if constexpr(containsForces()) {
                forces_[ix] = {};
            }
            if constexpr(containsVelocities()) {
                velocities_[ix] = {};
            }
            particleTypes_[ix] = type;
            blanks.pop_back();
        }
    }

    void removeParticle(size_type index) {
        positions_[index].reset();
    }

    template<typename IteratorAdd, typename IteratorRemove>
    void update(IteratorAdd beginAdd, IteratorAdd endAdd, IteratorRemove beginRemove, IteratorRemove endRemove) {
        auto itAdd = beginAdd;
        auto itRemove = beginRemove;
        while(itAdd != endAdd && itRemove != endRemove) {
            const auto &[pos, type] = *itAdd;
            positions_[*itRemove] = pos;
            particleTypes_[*itRemove] = type;
            if constexpr(containsForces()) {
                forces_[*itRemove] = {};
            }
            if constexpr(containsVelocities()) {
                velocities_[*itRemove] = {};
            }
            ++itAdd;
            ++itRemove;
        }
        while (itRemove != endRemove) {
            removeParticle(*itRemove);
            ++itRemove;
        }
        while (itAdd != endAdd) {
            addParticle(std::get<0>(*itAdd), std::get<1>(*itAdd));
            ++itAdd;
        }
    }

    [[nodiscard]] ParticleType typeOf(std::size_t ix) const {
        return particleTypes_[ix];
    }

    template<typename F, typename Pool,
            typename R=std::invoke_result_t<std::decay_t<F>, std::size_t, Position&, ParticleType, Force&>,
            typename = std::enable_if_t<std::is_void_v<R>>>
    std::vector<std::future<void>> forEachParticle(F &&op, config::PoolPtr<Pool> pool) {
        std::vector<std::future<void>> futures;
        auto granularity = config::threadGranularity(pool);
        futures.reserve(granularity);

        auto loop = [operation = std::forward<F>(op)](
                auto startIndex,
                const auto &beginPositions, const auto &endPositions,
                auto itTypes, auto itForces, auto itVelocities
        ) {
            for (auto itPos = beginPositions; itPos != endPositions; ++itPos, ++startIndex) {
                if (*itPos) {
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

    template<typename F, typename Pool,
            typename R=std::invoke_result_t<std::decay_t<F>, std::size_t, Position&, ParticleType, Force&>,
            typename = std::enable_if_t<!std::is_void_v<R>>>
    std::vector<std::future<std::vector<R>>> forEachParticle(F &&op, config::PoolPtr<Pool> pool) {
        std::vector<std::future<std::vector<R>>> futures;
        auto granularity = config::threadGranularity(pool);
        futures.reserve(granularity);

        auto loop = [operation = std::forward<F>(op)](
                auto startIndex,
                const auto &beginPositions, const auto &endPositions,
                auto itTypes, auto itForces, auto itVelocities
        ) {
            std::vector<R> result;
            result.reserve(std::distance(beginPositions, endPositions));
            for (auto itPos = beginPositions; itPos != endPositions; ++itPos, ++startIndex) {
                if (*itPos) {
                    if constexpr(containsForces() && containsVelocities()) {
                        result.push_back(operation(startIndex, **itPos, *itTypes, *itForces, *itVelocities));
                    } else if constexpr(containsForces() && !containsVelocities()) {
                        result.push_back(operation(startIndex, **itPos, *itTypes, *itForces));
                    } else if constexpr(!containsForces() && containsVelocities()) {
                        result.push_back(operation(startIndex, **itPos, *itTypes, *itVelocities));
                    } else {
                        result.push_back(operation(startIndex, **itPos, *itTypes));
                    }
                }

                if constexpr(containsForces()) {
                    ++itForces;
                }
                if constexpr(containsVelocities()) {
                    ++itVelocities;
                }
            }
            return result;
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

    const ContainerType<MaybePosition> &positions() const {
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
    std::vector<size_type> blanks;
};

template<typename ParticleCollection>
struct ParticleCollectionUpdater {
    using Particles = ParticleCollection;
    using ParticleType = typename ParticleCollection::ParticleType;
    using Position = typename ParticleCollection::Position;
    using Index = typename ParticleCollection::size_type;

    explicit ParticleCollectionUpdater(const ParticleCollection &collection) : changed(collection.size()) {}

    std::vector<std::atomic<bool>> changed;
    std::deque<std::tuple<Position, ParticleType>> toAdd;
    std::deque<Index> toRemove;

    std::mutex addMutex, removeMutex;

    void add(const ParticleType &type, const Position &position) {
        std::scoped_lock lock(addMutex);
        toAdd.push_back(std::make_tuple(position, type));
    }

    void remove(const Index &index, ParticleCollection &collection) {
        bool expected {false};
        if(changed[index].compare_exchange_weak(expected, true)) {
            std::scoped_lock lock(removeMutex);
            collection.removeParticle(index);
        }
    }

    void directUpdate(const Index &index, const std::optional<ParticleType> &type,
                      const std::optional<Position> &position,
                      ParticleCollection &collection) {
        bool expected {false};
        if(changed[index].compare_exchange_weak(expected, true)) {
            if(type) {
                collection.setType(index, *type);
            }
            if(position) {
                collection.setPosition(index, *position);
            }
        }
    }
};

}
