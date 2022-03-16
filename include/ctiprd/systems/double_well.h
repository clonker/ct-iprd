//
// Created by mho on 3/16/22.
//

#pragma once

#include <cstddef>
#include <tuple>

#include <ctiprd/potentials/external.h>

namespace ctiprd::systems {

template<typename T>
struct DoubleWell {

    using ExternalPotentials = std::tuple<potentials::DoubleWell>;
    using PairPotentials = std::tuple<>;



};
}