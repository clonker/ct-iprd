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
#include "ctiprd/util/pbc.h"
#include <ctiprd/cpu/ContainerContainer.h>

namespace ctiprd::cpu {

namespace detail {

template<typename Info, typename F>
struct Loop {

    explicit Loop(const F &f) : operation(std::cref(f)) {}

    template<typename... Args>
    void operator()(std::size_t startId, std::size_t n, Args&&... args) const {
        static_assert(sizeof...(args) == 1 + static_cast<int>(Info::positions) + static_cast<int>(Info::velocities) + static_cast<int>(Info::forces));
        std::tuple<Args...> iterators {std::forward<Args>(args)...};
        for(std::size_t id = startId; id < startId + n; ++id) {
            if(*std::get<0>(iterators)) {
                if constexpr(Info::forces && Info::velocities) {
                    auto &[itPos, itTypes, itForces, itVelocities] = iterators;
                    operation(id, **itPos, *itTypes, *itForces, *itVelocities);
                } else if constexpr(Info::forces && !Info::velocities) {
                    auto &[itPos, itTypes, itForces] = iterators;
                    operation(id, **itPos, *itTypes, *itForces);
                } else if constexpr(!Info::forces && Info::velocities) {
                    auto &[itPos, itTypes, itVelocities] = iterators;
                    operation(id, **itPos, *itTypes, *itVelocities);
                } else {
                    auto &[itPos, itTypes] = iterators;
                    operation(id, **itPos, *itTypes);
                }
            }

            /*std::apply([](auto&... it) {
                (++it, ...);
            }, iterators);*/
        }
    }

    std::reference_wrapper<const F> operation;

};
}

namespace particles {

struct forces {};
struct velocities {};
struct positions {};

template<typename... args>
struct Info {
    static constexpr bool forces = false;
    static constexpr bool velocities = false;
    static constexpr bool positions = false;
};

template<typename... rest> struct Info<forces, rest...> : Info<rest...> {
    static constexpr bool forces = true;
};

template<typename... rest> struct Info<velocities, rest...> : Info<rest...> {
    static constexpr bool velocities = true;
};

template<typename... rest> struct Info<positions, rest...> : Info<rest...> {
    static constexpr bool positions = true;
};


}

template<typename System, typename... Args>
class ParticleCollection {
public:
    using Info = particles::Info<Args...>;
    static_assert(Info::positions, "currently only implemented with positions.");

    using dtype = typename System::dtype;
    static constexpr std::size_t DIM = System::DIM;
    using SystemType = System;

    using Position = Vec<dtype, DIM>;
    using State = Position;
    using MaybePosition = std::optional<Position>;
    using Force = Vec<dtype, DIM>;
    using Velocity = Vec<dtype, DIM>;
    using ParticleType = std::size_t;

    template<typename T>
    using ContainerType = std::vector<T>;

    using size_type = typename ContainerType<MaybePosition>::size_type;

    static constexpr int dim = DIM;

    static constexpr bool containsPositions() { return Info::positions; }
    static constexpr bool containsForces() { return Info::forces; }
    static constexpr bool containsVelocities() { return Info::velocities; }

    [[nodiscard]] std::size_t nParticles() const {
        return positions_.size() - blanks.size();
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

    const Position &positionOf(size_type index) const {
        return *positions_[index];
    }

    template<typename T>
    void initializeParticles(std::size_t n, T &&type) requires std::convertible_to<T, std::string_view> {
        auto &generator = ctiprd::rnd::staticThreadLocalGenerator();
        for(std::size_t i = 0; i < n; ++i) {
            Position pos {};
            for (std::size_t d = 0; d < DIM; ++d) {
                std::uniform_real_distribution<float> dist {-System::boxSize[d] / 2, System::boxSize[d] / 2};
                pos[d] = dist(generator);
            }

            addParticle(pos, type);
        }
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
        blanks.push_back(index);
    }

    template<typename AddIterator, typename RemoveIterator>
    void update(AddIterator beginAdd, AddIterator endAdd, RemoveIterator beginRemove, RemoveIterator endRemove) {
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

    [[nodiscard]] ParticleType typeOf(const std::size_t particleIndex) const {
        return particleTypes_[particleIndex];
    }

    [[nodiscard]] bool exists(const std::size_t particleIndex) const {
        return positions_[particleIndex].has_value();
    }

    void sort() {
        std::sort(begin(blanks), end(blanks));
        std::size_t nSwapped {0};
        for(std::size_t i = size() - 1; i > size() - 1 - blanks.size(); --i) {
            if(exists(i)) {
                // swap with nSwapped-th blank
                auto blank = blanks[nSwapped];
                positions_[blank] = positions_[i];
                particleTypes_[blank] = particleTypes_[i];
                if constexpr(containsForces()) {
                    forces_[blank] = forces_[i];
                }
                if constexpr(containsVelocities()) {
                    velocities_[blank] = velocities_[i];
                }
                ++nSwapped;
            }
        }
        positions_.resize(positions_.size() - blanks.size());
        blanks.clear();
    }

    template<typename F, typename Pool>
    std::vector<std::future<void>> forEachParticle(F &&operation, config::PoolPtr<Pool> pool) {
        std::vector<std::future<void>> futures;
        const auto granularity = config::threadGranularity(pool);
        futures.reserve(granularity);

        const auto loop = [operation = std::forward<F>(operation)](
                auto startIndex,
                const auto &beginPositions, const auto &endPositions,
                auto itTypes, auto itForces, auto itVelocities
        ) {
            for (auto itPos = beginPositions; itPos != endPositions; ++itPos, ++startIndex, ++itTypes) {
                if (*itPos) {
                    if constexpr(containsForces() && containsVelocities()) {
                        operation(startIndex, **itPos, *itTypes, *itForces, *itVelocities);
                    } else if constexpr(containsForces() && !containsVelocities()) {
                        operation(startIndex, **itPos, *itTypes, *itForces);
                    } else if constexpr(!containsForces() && containsVelocities()) {
                        operation(startIndex, **itPos, *itTypes, +itVelocities);
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

        const auto grainSize = size() / granularity;
        auto itPos = positions_.begin();
        auto itTypes = particleTypes_.begin();
        auto itForces = forces_.begin();
        auto itVelocities = velocities_.begin();

        std::size_t beginIx = 0;
        for (std::size_t i = 0; i < granularity - 1; ++i) {
            if (itPos != positions_.end()) {
                futures.emplace_back(pool->push(loop, beginIx, itPos, itPos + grainSize, itTypes, itForces, itVelocities));
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
            futures.emplace_back(pool->push(loop, beginIx, itPos, end(positions_), itTypes, itForces, itVelocities));
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

template<typename T, typename ParticleCollection>
struct ParticleCollectionUpdater {
    using System = T;
    using Particles = ParticleCollection;
    using ParticleType = typename ParticleCollection::ParticleType;
    using Position = typename ParticleCollection::Position;
    using State = Position;
    using Index = typename ParticleCollection::size_type;
    using dtype = typename Particles::dtype;
    static constexpr int dim = Particles::dim;

    explicit ParticleCollectionUpdater(const ParticleCollection &collection) : changed(collection.size()) {}

    std::vector<std::atomic<bool>> changed;
    std::deque<std::tuple<Position, ParticleType>> toAdd;
    std::deque<Index> toRemove;

    // std::mutex addMutex, removeMutex;

    void add(const ParticleType &type, Position &&position, const auto &/*ignore*/) {
        util::pbc::wrapPBC<System>(position);
        // std::scoped_lock lock(addMutex);
        toAdd.emplace_back(std::move(position), type);
    }

    template<bool check=true>
    void remove(const Index &index, ParticleCollection &collection) {
        if constexpr(check) {
            bool expected{false};
            if (changed[index].compare_exchange_weak(expected, true)) {
                // std::scoped_lock lock(removeMutex);
                collection.removeParticle(index);
            }
        } else {
            // std::scoped_lock lock(removeMutex);
            collection.removeParticle(index);
        }
    }

    bool claim(const Index &index) {
        bool expected {false};
        return changed[index].compare_exchange_weak(expected, true);
    }

    void release(const Index &index) {
        changed[index] = false;
    }

    template<bool check=true>
    void directUpdate(const Index &index, const std::optional<ParticleType> &type,
                      std::optional<Position> &&position,
                      ParticleCollection &collection) {
        if constexpr(check) {
            bool expected {false};
            if(changed[index].compare_exchange_weak(expected, true)) {
                if(type) {
                    collection.setType(index, *type);
                }
                if(position) {
                    util::pbc::wrapPBC<System>(*position);
                    collection.setPosition(index, *position);
                }
            }
        } else {
            if(type) {
                collection.setType(index, *type);
            }
            if(position) {
                util::pbc::wrapPBC<System>(*position);
                collection.setPosition(index, *position);
            }
        }
    }
};

template<typename T, typename ParticleCollection>
struct SynchronizedParticleCollectionUpdater {
    using Particles = ParticleCollection;
    using ParticleType = typename ParticleCollection::ParticleType;
    using Position = typename ParticleCollection::Position;
    using State = Position;
    using Index = typename ParticleCollection::size_type;
    using dtype = typename Particles::dtype;
    using System = T;
    static constexpr int dim = Particles::dim;
    static constexpr bool periodic = T::periodic;

    explicit SynchronizedParticleCollectionUpdater(const ParticleCollection &collection) : changed(collection.size()) {}

    std::vector<char> changed;

    void add(const ParticleType &type, Position &&position, ParticleCollection &collection) {
        util::pbc::wrapPBC<System>(position);
        collection.addParticle(position, type);
    }

    template<bool check=true>
    void remove(const Index &index, ParticleCollection &collection) {
        if constexpr(check) {
            if(!changed[index]) {
                changed[index] = true;
                collection.removeParticle(index);
            }
        } else {
            collection.removeParticle(index);
        }
    }

    bool claim(const Index &index) {
        if(!changed[index]) {
            changed[index] = true;
            return true;
        }
        return false;
    }

    void release(const Index &index) {
        changed[index] = false;
    }

    template<bool check=true>
    void directUpdate(const Index &index, const std::optional<ParticleType> &type,
                      std::optional<Position> &&position,
                      ParticleCollection &collection) {
        if constexpr(check) {
            bool expected {false};
            if(!changed[index]) {
                changed[index] = true;
                if(type) {
                    collection.setType(index, *type);
                }
                if(position) {
                    util::pbc::wrapPBC<System>(*position);
                    collection.setPosition(index, *position);
                }
            }
        } else {
            if(type) {
                collection.setType(index, *type);
            }
            if(position) {
                util::pbc::wrapPBC<System>(*position);
                collection.setPosition(index, *position);
            }
        }
    }
};

}
