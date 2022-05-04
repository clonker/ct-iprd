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
#include <ctiprd/util/rates.h>

namespace ctiprd::systems {

template<typename T>
struct Conf {
    static constexpr T diffPrey = 0.01;
    static constexpr T diffPred = 0.01;

    static constexpr T alpha = 2.;  // birth: prey -> prey + prey
    static constexpr T alphaDistance = 2.;  // birth distance
    static constexpr T beta = 0.05;  // eat: prey + pred -> pred + pred
    static constexpr T betaRadius = 0.25;
    static constexpr T betaMic = 7.670679846561291;
    static constexpr T gamma = 1.5;

    static constexpr T friction = 0.01;
    static constexpr T frictionRadius = 0.2;
    static constexpr T frictionMic = 0.39155565024786165;
};

template<typename T>
struct LotkaVolterra {
    using Cfg = Conf<T>;
    using dtype = T;
    static constexpr std::size_t DIM = 3;
    static constexpr std::array<T, DIM> boxSize{5., 5., 5.};
    static constexpr bool periodic = true;
    static constexpr T kBT = 2.43614;
    static constexpr ParticleTypes<dtype, 2> types{{
              {
                      .name = "predator",
                      .diffusionConstant = Cfg::diffPred
              },
              {
                      .name = "prey",
                      .diffusionConstant = Cfg::diffPrey
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
                    .productDistance = Cfg::alphaDistance,
                    .rate = Cfg::alpha
            };
            death = reactions::doi::Decay<T>{
                    .eductType = predatorId,
                    .rate = Cfg::gamma
            };
        }
        {
            auto &[preySocialFriction, predatorSocialFriction, predEatsPrey] = reactionsO2;
            preySocialFriction = reactions::doi::Fusion<T>{
                    .eductType1 = preyId,
                    .eductType2 = preyId,
                    .productType = preyId,
                    .reactionRadius = Cfg::frictionRadius,
                    .rate = Cfg::frictionMic
            };
            predatorSocialFriction = reactions::doi::Fusion<T>{
                    .eductType1 = predatorId,
                    .eductType2 = predatorId,
                    .productType = predatorId,
                    .reactionRadius = Cfg::frictionRadius,
                    .rate = Cfg::frictionMic
            };
            predEatsPrey = reactions::doi::Catalysis<T>{
                    .catalyst = predatorId,
                    .eductType = preyId,
                    .productType = predatorId,
                    .reactionRadius = Cfg::betaRadius,
                    .rate = Cfg::betaMic
            };
        }


        {
            static constexpr dtype eps = 1e-5;

            auto kmac = rates::macroscopicRate(Cfg::betaMic, types[preyId].diffusionConstant,
                                               types[predatorId].diffusionConstant, Cfg::betaRadius);
            if(std::abs(kmac - Cfg::beta) > eps) {
                throw std::runtime_error(fmt::format("kmac = {}, beta = {}", kmac, Cfg::beta));
            }
        }
        {
            static constexpr dtype eps = 1e-5;

            auto kmac = rates::macroscopicRate(Cfg::frictionMic, types[preyId].diffusionConstant,
                                               types[predatorId].diffusionConstant, Cfg::frictionRadius);
            if(std::abs(kmac - Cfg::friction) > eps) {
                throw std::runtime_error(fmt::format("kmac = {}, friction = {}", kmac, Cfg::friction));
            }
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
