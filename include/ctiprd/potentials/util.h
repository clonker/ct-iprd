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

}
