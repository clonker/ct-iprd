//
// Created by mho on 3/16/22.
//
#pragma once

#include <cstddef>
#include <cmath>

namespace ctiprd::potentials::external {

template<typename dtype, std::size_t typeId>
struct DoubleWell {

    constexpr static std::size_t particleType = typeId;

    static constexpr std::size_t DIM = 2;

    template<typename State, typename ParticleType>
    [[nodiscard]] constexpr auto energy(const State &x, const ParticleType &type) const {
        if (type == particleType) {
            return k*(x[0] * x[0] - 1.) * (x[0] * x[0] - 1.) + k*x[1] * x[1];
        } else {
            return static_cast<dtype>(0);
        }
    }

    template<typename State, typename ParticleType>
    [[nodiscard]] constexpr State force(const State &x, const ParticleType &type) const {
        if (type == particleType) {
            return {{-4 * k * x[0] * x[0] * x[0] + 4 * k * x[0], -2 * k * x[1]}};
        } else {
            return {};
        }
    }

    dtype k {1.};
};

}
