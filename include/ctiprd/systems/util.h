//
// Created by mho on 3/25/22.
//

#pragma once

#include <string_view>
#include <algorithm>

namespace ctiprd::systems {

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

template<typename System>
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

// todo would be better to do this for System!
template<typename T>
concept system_info = requires {
    { T::SystemType };
    { T::dtype };
    { T::periodic } -> std::same_as<bool>;
    { T::DIM } -> std::same_as<std::size_t>;
    { T::boxSize };  // todo check for array and dims
    { T::nTypes } -> std::same_as<std::size_t>;
    { T::ExternalPotentials };
    { T::PairPotentials };
    { T::hasForces } -> std::same_as<bool>;
    // todo check rest
};

}
