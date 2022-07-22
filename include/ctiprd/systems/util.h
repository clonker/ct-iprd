//
// Created by mho on 3/25/22.
//

#pragma once

#include <string_view>
#include <algorithm>

namespace ctiprd::systems {

template<typename T>
concept system = requires {
    { T::dtype };
    /*{ T::periodic };
    { T::DIM };
    { T::kBT };
    { T::boxSize };  // todo check for array and dims
    { T::ExternalPotentials };
    { T::PairPotentials };
    { T::ReactionsO1 };
    { T::ReactionsO2 };
    { T::types };*/
};

template<auto &types>
static constexpr std::size_t particleTypeId(std::string_view name) {
    std::size_t i = 0;
    for(auto it = begin(types); it != end(types); ++it, ++i) {
        if (it->name == name) {
            break;
        }
    }
    return i;
}

template<system System>
struct SystemInfo {
    using SystemType = System;
    using dtype = typename System::dtype;
    static constexpr bool periodic = System::periodic;
    static constexpr std::size_t DIM = System::DIM;
    static constexpr std::array<dtype, DIM> boxSize = System::boxSize;
    static constexpr std::size_t nTypes = System::types.size();

    using ExternalPotentials = typename System::ExternalPotentials;
    using PairPotentials = typename System::PairPotentials;

    static constexpr int nExternalPotentials = std::tuple_size_v<ExternalPotentials>;
    static constexpr int nPairPotentials = std::tuple_size_v<PairPotentials>;

    static constexpr bool hasForces() {
        return nExternalPotentials > 0 || nPairPotentials > 0;
    }

    using ReactionsO1 = typename System::ReactionsO1;
    using ReactionsO2 = typename System::ReactionsO2;

    static constexpr int nReactionsO1 = std::tuple_size_v<ReactionsO1>;
    static constexpr int nReactionsO2 = std::tuple_size_v<ReactionsO2>;

    static constexpr bool hasReactions() {
        return nReactionsO1 > 0 || nReactionsO2 > 0;
    }

    static constexpr std::size_t particleId(std::string_view name) {
        return particleTypeId<System::types>(name);
    }

    static constexpr dtype diffusionConstantOf(std::size_t typeId) {
        return System::types[typeId].diffusionConstant;
    }

    static constexpr dtype diffusionConstantOf(std::string_view name) {
        return diffusionConstantOf(particleId(name));
    }
};

}
