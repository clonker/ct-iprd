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
#include <ctiprd/potentials/interaction.h>
#include <ctiprd/ParticleTypes.h>
#include <ctiprd/systems/util.h>

namespace ctiprd::systems {

template<typename T>
struct DoubleWell {
    using dtype = T;
    static constexpr std::size_t DIM = 2;
    static constexpr std::array<T, DIM> boxSize {5., 5.};

    static constexpr ParticleTypes<dtype, 3> types {{
            {
                .name = "A",
                .diffusionConstant = 1.
            },
            {
                .name = "B",
                .diffusionConstant = 1.
            },
            {
                .name = "C",
                .diffusionConstant = 55.
            }
    }};

    using ExternalPotentials = std::tuple<
            potential::external::DoubleWell<T, particleTypeId<types>("B")>,
            potential::external::DoubleWell<T, particleTypeId<types>("C")>
    >;
    using PairPotentials = std::tuple<
            potential::pair::HarmonicRepulsion<T, particleTypeId<types>("A"), particleTypeId<types>("B")>
    >;
};
}
