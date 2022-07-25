//
// Created by mho on 3/25/22.
//

#pragma once

#include <string_view>
#include <algorithm>
#include <type_traits>

#include <ctiprd/ParticleTypes.h>

namespace ctiprd::systems {

template<typename T>
concept system = requires (T instance){
    { typename T::dtype{} } -> std::floating_point;
    { T::periodic } -> std::convertible_to<bool>;
    { T::DIM } -> std::convertible_to<std::size_t>;
    { T::kBT } -> std::convertible_to<typename T::dtype>;
    { T::boxSize } -> std::convertible_to<std::array<typename T::dtype, T::DIM>>;
    { instance.externalPotentials };
    { instance.pairPotentials };
    { instance.reactionsO1 };
    { instance.reactionsO2 };
    { instance.types };
    { std::tuple_size_v<typename T::ExternalPotentials> };
    { std::tuple_size_v<typename T::PairPotentials> };
    { std::tuple_size_v<typename T::ReactionsO1> };
    { std::tuple_size_v<typename T::ReactionsO2> };
    { instance.types[0] } -> std::convertible_to<ParticleType<typename T::dtype>>;
};

template<auto &types>
static constexpr std::size_t particleTypeId(std::string_view name) {
    std::size_t typeIndex = 0;
    for(const auto &type : types) {
        if (type.name == name) { break; }
        ++typeIndex;
    }
    return typeIndex;
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
