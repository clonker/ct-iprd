//
// Created by mho on 3/16/22.
//

#pragma once

#include <cstddef>
#include <tuple>
#include <optional>

#include <ctiprd/vec.h>
#include <ctiprd/integrators/EulerMaruyama.h>
#include <ctiprd/potentials/external.h>

namespace ctiprd::systems {

template<typename T>
struct DoubleWell {
    using State = Vec<T, 2>;

    using ExternalPotentials = std::tuple<
            potential::external::DoubleWell<State, State>
    >;
    using PairPotentials = std::tuple<>;

    using Integrator = integrator::EulerMaruyama<State::dim, T, ExternalPotentials, PairPotentials>;
};
}
