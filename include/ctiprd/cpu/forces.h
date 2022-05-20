//
// Created by mho on 5/20/22.
//

#pragma once

#include <memory>
#include <list>
#include <vector>

#include <tsl/robin_map.h>

#include "ctiprd/util/hash.h"
#include "ctiprd/potentials/external.h"
#include "ctiprd/potentials/interaction.h"
#include "ctiprd/systems/util.h"
#include "ctiprd/util/ops.h"

namespace ctiprd::cpu::forces {

template<typename ParticleCollection>
struct ForceO1 {
    using dtype = typename ParticleCollection::dtype;
    using State = typename ParticleCollection::State;

    [[nodiscard]] virtual dtype energy([[maybe_unused]] const State &state) const { return 0; }
    [[nodiscard]] virtual State force([[maybe_unused]] const State &state) const { return {}; }
};
template<typename ParticleCollection>
struct ForceO2 {
    using dtype = typename ParticleCollection::dtype;
    using State = typename ParticleCollection::State;

    [[nodiscard]] virtual dtype energy([[maybe_unused]] const State &state1,
                                       [[maybe_unused]] const State &state2) const { return 0; }
    [[nodiscard]] virtual State force([[maybe_unused]] const State &state1,
                                      [[maybe_unused]] const State &state2) const { return {}; }
};

namespace detail {

template<typename ParticleCollection, typename Potential>
struct CPUForceO1 : public ForceO1<ParticleCollection> {};

template<typename ParticleCollection, typename Potential>
struct CPUForceO2 : public ForceO2<ParticleCollection> {};

namespace geom {

template<bool inclusion, typename Geometry, typename State, typename dtype>
dtype energy(const dtype &forceConstant, const State &state, const Geometry &geometry) {
    auto shortestDiff = geometry.template smallestDifference<inclusion>(state);
    return forceConstant * shortestDiff.normSquared() / 2;
}

template<bool inclusion, typename Geometry, typename State, typename dtype>
State force(const dtype &forceConstant, const State &state, const Geometry &geometry) {
    return -1 * forceConstant * geometry.template smallestDifference<inclusion>(state);
}

}

template<typename ParticleCollection>
struct CPUForceO1<ParticleCollection, ctiprd::potentials::external::BoxInclusion<typename ParticleCollection::dtype, ParticleCollection::DIM>> : ForceO1<ParticleCollection> {
    using typename ForceO1<ParticleCollection>::dtype;
    using typename ForceO1<ParticleCollection>::State;
    using BasePotential = ctiprd::potentials::external::BoxInclusion<typename ParticleCollection::dtype, ParticleCollection::DIM>;

    explicit CPUForceO1(const BasePotential &basePotential) : basePotential(basePotential) {}

    [[nodiscard]] dtype energy(const State &state) const override {
        return geom::energy<true>(basePotential.k, state, basePotential.geometry);
    }

    [[nodiscard]] State force(const State &state) const override {
        return geom::force<true>(basePotential.k, state, basePotential.geometry);
    }

    BasePotential basePotential;
};

template<typename ParticleCollection>
struct CPUForceO1<ParticleCollection, ctiprd::potentials::external::BoxExclusion<typename ParticleCollection::dtype, ParticleCollection::DIM>> : ForceO1<ParticleCollection> {
    using typename ForceO1<ParticleCollection>::dtype;
    using typename ForceO1<ParticleCollection>::State;
    using BasePotential = ctiprd::potentials::external::BoxExclusion<typename ParticleCollection::dtype, ParticleCollection::DIM>;

    explicit CPUForceO1(const BasePotential &basePotential) : basePotential(basePotential) {}

    [[nodiscard]] dtype energy(const State &state) const override {
        return geom::energy<false>(basePotential.k, state, basePotential.geometry);
    }

    [[nodiscard]] State force(const State &state) const override {
        return geom::force<false>(basePotential.k, state, basePotential.geometry);
    }

    BasePotential basePotential;
};

template<typename ParticleCollection>
struct CPUForceO1<ParticleCollection, ctiprd::potentials::external::DoubleWell<typename ParticleCollection::dtype>> : ForceO1<ParticleCollection> {
    using typename ForceO1<ParticleCollection>::dtype;
    using typename ForceO1<ParticleCollection>::State;
    using BasePotential = ctiprd::potentials::external::DoubleWell<typename ParticleCollection::dtype>;

    explicit CPUForceO1(const BasePotential &basePotential) : basePotential(basePotential) {}

    dtype energy(const State &state) const override {
        return basePotential.k*(state[0] * state[0] - 1.) * (state[0] * state[0] - 1.) + basePotential.k*state[1] * state[1];
    }

    State force(const State &state) const override {
        return {{-4 * basePotential.k * state[0] * state[0] * state[0] + 4 * basePotential.k * state[0],
                 -2 * basePotential.k * state[1]}};
    }

    BasePotential basePotential;
};

template<typename ParticleCollection>
struct CPUForceO2<ParticleCollection, ctiprd::potentials::pair::HarmonicRepulsion<typename ParticleCollection::dtype>> : ForceO2<ParticleCollection> {
    using typename ForceO2<ParticleCollection>::dtype;
    using typename ForceO2<ParticleCollection>::State;
    using BasePotential = ctiprd::potentials::pair::HarmonicRepulsion<typename ParticleCollection::dtype>;

    explicit CPUForceO2(const BasePotential &basePotential) : basePotential(basePotential),
                    cutoffSquared(basePotential.cutoff * basePotential.cutoff) {}

    dtype energy(const State &state1, const State &state2) const override {
        auto dSquared = ctiprd::util::pbc::dSquared<typename ParticleCollection::SystemType>(state1, state2);
        if (dSquared < cutoffSquared) {
            auto dist = std::sqrt(dSquared);
            dist -= basePotential.cutoff;
            return dist * dist * basePotential.forceConstant / 2;
        }
        return 0;
    }

    State force(const State &state1, const State &state2) const override {
        auto xij = util::pbc::shortestDifference<typename ParticleCollection::SystemType>(state1, state2);
        auto dSquared = xij.normSquared();
        if (dSquared < cutoffSquared && dSquared > 0) {
            auto dist = std::sqrt(dSquared);
            return (basePotential.forceConstant * (dist - basePotential.cutoff)) / dist * xij;
        }

        return {};
    }

private:
    BasePotential basePotential;
    dtype cutoffSquared;
};

using KeyO2 = std::tuple<std::size_t, std::size_t>;
using HasherO2 = ctiprd::util::hash::ForwardBackwardTupleHasher<KeyO2>;
using EqO2 = ctiprd::util::hash::ForwardBackwardTupleEquality<KeyO2>;
}

template<typename ParticleCollection>
using FFO1Map = tsl::robin_map<std::size_t, std::vector<const ForceO1<ParticleCollection>*>>;
template<typename ParticleCollection>
using FFO1Backing = std::list<std::unique_ptr<ForceO1<ParticleCollection>>>;

template<typename ParticleCollection>
using FFO2Map = tsl::robin_map<detail::KeyO2, std::vector<const ForceO2<ParticleCollection>*>, detail::HasherO2, detail::EqO2>;
template<typename ParticleCollection>
using FFO2Backing = std::list<std::unique_ptr<ForceO2<ParticleCollection>>>;

template<typename ParticleCollection>
std::tuple<FFO1Map<ParticleCollection>, FFO1Backing<ParticleCollection>> generateMapO1(const typename ParticleCollection::SystemType &system) {
    using System = typename ParticleCollection::SystemType;
    using SystemInfo = systems::SystemInfo<System>;
    FFO1Backing<ParticleCollection> backingData {};
    FFO1Map<ParticleCollection> map;
    map.reserve(SystemInfo::nExternalPotentials);

    // allocate for each type
    for (const auto &[name, _] : System::types) {
        map[ctiprd::systems::particleTypeId<System::types>(name)] = {};
    }

    // now push back implementations
    [&]<auto... I>(std::index_sequence<I...>) {
        ([&](const auto &potential) {
            using ForceType = std::tuple_element_t<I, typename System::ExternalPotentials>;
            backingData.push_back(std::make_unique<detail::CPUForceO1<ParticleCollection, ForceType>>(potential));
            const auto &ref = backingData.back();

            for (const auto &[name, _] : System::types) {
                if(potential.supportsType(ctiprd::systems::particleTypeId<System::types>(name))) {
                    map[ctiprd::systems::particleTypeId<System::types>(name)].push_back(ref.get());
                }
            }
        }(std::get<I>(system.externalPotentials)), ...);
    }(std::make_index_sequence<std::tuple_size_v<typename System::ExternalPotentials>>{});
    return std::make_tuple(map, std::move(backingData));
}

template<typename ParticleCollection>
std::tuple<FFO2Map<ParticleCollection>, FFO2Backing<ParticleCollection>> generateMapO2(const typename ParticleCollection::SystemType &system) {
    using System = typename ParticleCollection::SystemType;
    using SystemInfo = typename systems::SystemInfo<System>;
    FFO2Backing<ParticleCollection> backingData {};
    FFO2Map<ParticleCollection> map;
    map.reserve(ctiprd::util::binomialCoefficient(SystemInfo::nPairPotentials, 2));

    {
        for (const auto &[name1, d1] : System::types) {
            auto id1 = ctiprd::systems::particleTypeId<System::types>(name1);
            for (const auto &[name2, d2] : System::types) {
                auto id2 = ctiprd::systems::particleTypeId<System::types>(name2);
                map[{id1, id2}] = {};
            }
        }
    }

    [&]<auto... I>(std::index_sequence<I...>) {
        ([&](const auto &potential) {
            using PotentialType = std::tuple_element_t<I, typename System::PairPotentials>;
            using PotentialImpl = detail::CPUForceO2<ParticleCollection, PotentialType>;
            backingData.push_back(std::make_unique<PotentialImpl>(potential));
            const auto &ref = backingData.back();

            for (const auto &[name1, d1] : System::types) {
                auto id1 = ctiprd::systems::particleTypeId<System::types>(name1);
                for (const auto &[name2, d2] : System::types) {
                    auto id2 = ctiprd::systems::particleTypeId<System::types>(name2);
                    if(potential.supportsTypes(id1, id2)) {
                        map[std::make_tuple(id1, id2)].push_back(ref.get());
                    }
                }
            }
        }(std::get<I>(system.pairPotentials)), ...);
    }(std::make_index_sequence<std::tuple_size_v<typename System::PairPotentials>>{});
    return std::make_tuple(map, std::move(backingData));
}

}
