//
// Created by mho on 4/5/22.
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
struct LotkaVolterra {
    using dtype = T;
    static constexpr std::size_t DIM = 2;
    static constexpr std::array<T, DIM> boxSize{10., 50.};
    static constexpr bool periodic = false;
    static constexpr ParticleTypes<dtype, 2> types{{
              {
                      .name = "predator",
                      .diffusionConstant = .1
              },
              {
                      .name = "prey",
                      .diffusionConstant = .1
              },
    }};
    static constexpr std::size_t preyId = particleTypeId<types>("prey");
    static constexpr std::size_t predatorId = particleTypeId<types>("predator");
    
    LotkaVolterra() {
        {
            std::get<0>(reactionsO1) = reactions::doi::Fission<T>{
                    .eductType = preyId,
                    .productType1 = preyId,
                    .productType2 = preyId,
                    .productDistance = .5,
                    .rate = 2
            };
            std::get<1>(reactionsO1) = reactions::doi::Decay<T>{
                    .eductType = predatorId,
                    .rate = 1.5
            };
        }
        {
            std::get<0>(reactionsO2) = reactions::doi::Fusion<T>{
                    .eductType1 = preyId,
                    .eductType2 = preyId,
                    .productType = preyId,
                    .reactionRadius = 0.05,
                    .rate = 21.112150808892327
            };
            std::get<1>(reactionsO2) = reactions::doi::Fusion<T>{
                    .eductType1 = predatorId,
                    .eductType2 = predatorId,
                    .productType = predatorId,
                    .reactionRadius = 0.05,
                    .rate = 21.112150808892327
            };
            std::get<2>(reactionsO2) = reactions::doi::Catalysis<T>{
                    .catalyst = predatorId,
                    .eductType = preyId,
                    .productType = predatorId,
                    .reactionRadius = 0.3,
                    .rate = 19.371705844861612
            };
        }
        {
            auto &box = std::get<0>(externalPotentials);
            box.geometry.v0 = {-4., -24.};
            box.geometry.v1 = {4., 24.};
            box.k = 50.;
        }
    }


    using ExternalPotentials = std::tuple<potentials::external::BoxInclusion<dtype, DIM, true>>;
    using PairPotentials = std::tuple<>;

    using ReactionsO1 = std::tuple<
            reactions::doi::Fission<T>, // prey birth
            reactions::doi::Decay<T> // predator death
    >;
    using ReactionsO2 = std::tuple<
            reactions::doi::Fusion<T>, // prey social friction
            reactions::doi::Fusion<T>, // predator social friction
            reactions::doi::Catalysis<T> // predator eats prey
    >;

    ReactionsO1 reactionsO1{};
    ReactionsO2 reactionsO2{};

    ExternalPotentials externalPotentials{};
    PairPotentials pairPotentials{};
};

}
