//
// Created by mho on 3/16/22.
//
#pragma once

#include <cstddef>
#include <cmath>

namespace ctiprd::potentials::external {

template<typename dtype, bool allTypes, std::size_t typeId = 0>
struct DoubleWell {

    constexpr static std::size_t particleType = typeId;

    static constexpr std::size_t DIM = 2;

    template<typename State, typename ParticleType>
    [[nodiscard]] constexpr dtype energy(const State &x, const ParticleType &type) const {
        if (allTypes || type == particleType) {
            return k*(x[0] * x[0] - 1.) * (x[0] * x[0] - 1.) + k*x[1] * x[1];
        }
        return {};
    }

    template<typename State, typename ParticleType>
    [[nodiscard]] constexpr State force(const State &x, const ParticleType &type) const {
        if (allTypes || type == particleType) {
            return {{-4 * k * x[0] * x[0] * x[0] + 4 * k * x[0], -2 * k * x[1]}};
        }
        return {};
    }

    dtype k {1.};
};

}
