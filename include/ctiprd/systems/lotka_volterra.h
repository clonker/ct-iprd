//
// Created by mho on 4/5/22.
//

#pragma once

#include <cstddef>
#include <tuple>
#include <optional>

#include <ctiprd/vec.h>
#include <ctiprd/potentials/external.h>
#include <ctiprd/potentials/interaction.h>
#include <ctiprd/ParticleTypes.h>
#include <ctiprd/reactions/doi.h>
#include <ctiprd/systems/util.h>

namespace ctiprd::systems {

template<typename T>
struct LotkaVolterra {
    using dtype = T;
    static constexpr std::size_t DIM = 2;
    static constexpr std::array<T, DIM> boxSize{10., 10.};
    static constexpr bool periodic = true;
    static constexpr T kBT = 2.43614;
    static constexpr ParticleTypes<dtype, 2> types{{
              {
                      .name = "predator",
                      .diffusionConstant = 10.
              },
              {
                      .name = "prey",
                      .diffusionConstant = 10.
              },
    }};
    static constexpr std::size_t preyId = particleTypeId<types>("prey");
    static constexpr std::size_t predatorId = particleTypeId<types>("predator");
    
    LotkaVolterra() {
        {
            auto& [birth, death] = reactionsO1;
            birth = reactions::doi::Fission<T>{
                    .eductType = preyId,
                    .productType1 = preyId,
                    .productType2 = preyId,
                    .productDistance = 1.,
                    .rate = 2
            };
            death = reactions::doi::Decay<T>{
                    .eductType = predatorId,
                    .rate = 1.5
            };
        }
        {
            auto &[preySocialFriction, predatorSocialFriction, predEatsPrey] = reactionsO2;
            preySocialFriction = reactions::doi::Fusion<T>{
                    .eductType1 = preyId,
                    .eductType2 = preyId,
                    .productType = preyId,
                    .reactionRadius = 0.5,
                    .rate = 0.019116848082565446
            };
            predatorSocialFriction = reactions::doi::Fusion<T>{
                    .eductType1 = predatorId,
                    .eductType2 = predatorId,
                    .productType = predatorId,
                    .reactionRadius = 0.5,
                    .rate = 0.019116848082565446
            };
            predEatsPrey = reactions::doi::Catalysis<T>{
                    .catalyst = predatorId,
                    .eductType = preyId,
                    .productType = predatorId,
                    .reactionRadius = 0.5,
                    .rate = 0.09595107151845629
            };
        }
        {/*
            auto &box = std::get<0>(externalPotentials);
            box.geometry.v0 = {-4., -24.};
            box.geometry.v1 = {4., 24.};
            box.k = 50.;
        */}
    }


    using ExternalPotentials = std::tuple<>;  // potentials::external::BoxInclusion<dtype, DIM, true>
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
