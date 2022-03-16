//
// Created by mho on 3/16/22.
//
#pragma once

#include <cstddef>

namespace ctiprd::potential::external {

template<typename Position, typename Force>
struct DoubleWell {

    static constexpr std::size_t DIM = 2;

    [[nodiscard]] constexpr auto energy(const Position &x) const {
        return (x[0] * x[0] - 1.) * (x[0] * x[0] - 1.) + x[1] * x[1];
    }

    [[nodiscard]] constexpr Force force(const Position &x) const {
        return {{-4 * x[0] * x[0] * x[0] + 4 * x[0], -2 * x[1]}};
    }
};

}
