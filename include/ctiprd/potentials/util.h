//
// Created by mho on 3/25/22.
//
#pragma once

#include <tuple>

namespace ctiprd::potentials {

template<typename dtype, typename PairPotentials>
dtype cutoff(const PairPotentials &potentials) {
    dtype result {};
    std::apply([&result](auto &&... args) { ((result = std::max(result, args.cutoff)), ...); }, potentials);
    return result;
}

template<typename Potential, typename ParticleType>
bool potentialApplicable(const ParticleType &type) {
    return Potential::particleType == type;
}

template<typename Potential, typename ParticleType>
bool potentialApplicable(const ParticleType &t1, const ParticleType &t2) {
    return (Potential::particleType1 == t1 && Potential::particleType2 == t2)
        || (Potential::particleType1 == t2 && Potential::particleType2 == t1);
}

}
