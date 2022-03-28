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
#include <ctiprd/reactions/basic.h>
#include <ctiprd/systems/util.h>

namespace ctiprd::systems {

template<typename T>
struct DoubleWell {
    using dtype = T;
    static constexpr std::size_t DIM = 2;
    static constexpr std::array<T, DIM> boxSize {5., 5.}; // todo this should probably go into the integrator? or elsewhere? (for working w/ thermostat)
    static constexpr bool periodic = true;

    static constexpr ParticleTypes<dtype, 3> types {{
            {
                .name = "A",
                .diffusionConstant = 1.
            },
            {
                .name = "B",
                .diffusionConstant = 1.
            },
    }};

    using ExternalPotentials = std::tuple<
            potentials::external::DoubleWell<T, false, particleTypeId<types>("A")>,
            potentials::external::DoubleWell<T, false, particleTypeId<types>("B")>
    >;
    using PairPotentials = std::tuple<
            potentials::pair::HarmonicRepulsion<T, true>
    >;

    using ReactionsO1 = std::tuple<
            reactions::doi::Conversion<T>,
            reactions::doi::Conversion<T>
    >;
    using ReactionsO2 = std::tuple<
            reactions::doi::Catalysis<T>
    >;

    ReactionsO1 reactionsO1 {
            {
                .eductType = particleTypeId<types>("A"),
                .productType = particleTypeId<types>("B"),
                .rate = 1.
            },
            {
                .eductType = particleTypeId<types>("B"),
                .productType = particleTypeId<types>("B"),
                .rate = 1.
            }
    };
    ReactionsO2 reactionsO2 {
            {
                .rate = 1.,
                .eductRadius = .1,
            },
    };
};
}
