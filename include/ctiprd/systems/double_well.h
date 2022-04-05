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

    DoubleWell() {
        {
            auto &[conv1, conv2, decay] = reactionsO1;
            conv1.eductType = particleTypeId<types>("A");
            conv1.productType = particleTypeId<types>("B");
            conv1.rate = .0;

            conv2.eductType = particleTypeId<types>("B");
            conv2.productType = particleTypeId<types>("A");
            conv2.rate = 0;

            decay.eductType = particleTypeId<types>("A");
            decay.rate = 0;
        }
        {
            auto &[cata, fusion] = reactionsO2;
            cata.eductType = particleTypeId<types>("A");
            cata.catalyst = particleTypeId<types>("A");
            cata.productType = particleTypeId<types>("B");
            cata.rate = 0;
            cata.reactionRadius = .1;

            fusion.eductType1 = particleTypeId<types>("A");
            fusion.eductType2 = particleTypeId<types>("A");
            fusion.productType = particleTypeId<types>("B");
            fusion.rate = .1;
            fusion.reactionRadius = .1;
        }
    }

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

    /*using ReactionsO1 = std::tuple<
            reactions::doi::Conversion<T>,
            reactions::doi::Conversion<T>
    >;*/
    using ReactionsO1 = std::tuple<
            reactions::doi::Conversion<T>,
            reactions::doi::Conversion<T>,
            reactions::doi::Decay<T>
    >;
    using ReactionsO2 = std::tuple<
            reactions::doi::Catalysis<T>,
            reactions::doi::Fusion<T>
    >;

    ReactionsO1 reactionsO1 {};
    ReactionsO2 reactionsO2 {};
};
}
