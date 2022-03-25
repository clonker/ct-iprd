//
// Created by mho on 3/16/22.
//
#pragma once

#include <cstddef>
#include <cmath>

namespace ctiprd::potential::external {

template<typename dtype, std::size_t typeId>
struct DoubleWell {

    constexpr static std::size_t particleType = typeId;

    static constexpr std::size_t DIM = 2;
    using State = Vec<dtype, DIM>;

    [[nodiscard]] constexpr auto energy(const State &x) const {
        return k*(x[0] * x[0] - 1.) * (x[0] * x[0] - 1.) + k*x[1] * x[1];
    }

    [[nodiscard]] constexpr State force(const State &x) const {
        return {{-4 * k * x[0] * x[0] * x[0] + 4 * k * x[0], -2 * k * x[1]}};
    }

    dtype k {1.};
};

}
