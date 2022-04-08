/**
 *
 *
 * @file michaelis_menten.h
 * @brief 
 * @author clonker
 * @date 4/8/22
 */
#pragma once

#include <tuple>

#include <ctiprd/reactions/basic.h>
#include <ctiprd/systems/util.h>
#include <ctiprd/ParticleTypes.h>

namespace ctiprd::systems {
template<typename T>
struct MichaelisMenten {
    using dtype = T;
    static constexpr std::size_t DIM = 3;
    static constexpr std::array<T, DIM> boxSize{.3, .3, .3};
    static constexpr bool periodic = true;
    static constexpr T kBT = 2.43614;
    static constexpr ctiprd::ParticleTypes<dtype, 4> types{{
                                                                   { .name = "E", .diffusionConstant = 10.},
                                                                   { .name = "S", .diffusionConstant = 10.},
                                                                   { .name = "ES", .diffusionConstant = 10.},
                                                                   { .name = "P", .diffusionConstant = 10.}
                                                           }};
    static constexpr std::size_t eId = ctiprd::systems::particleTypeId<types>("E");
    static constexpr std::size_t sId = ctiprd::systems::particleTypeId<types>("S");
    static constexpr std::size_t esId = ctiprd::systems::particleTypeId<types>("ES");
    static constexpr std::size_t pId = ctiprd::systems::particleTypeId<types>("P");

    MichaelisMenten() {
        {
            std::get<0>(reactionsO1) = ctiprd::reactions::doi::Fission<T>{
                    .eductType = esId,
                    .productType1 = eId,
                    .productType2 = sId,
                    .productDistance = .03,
                    .rate = 1.
            };
            std::get<0>(reactionsO1) = ctiprd::reactions::doi::Fission<T>{
                    .eductType = esId,
                    .productType1 = eId,
                    .productType2 = pId,
                    .productDistance = .03,
                    .rate = 1.
            };
        }
        {
            std::get<0>(reactionsO2) = ctiprd::reactions::doi::Fusion<T>{
                    .eductType1 = eId,
                    .eductType2 = sId,
                    .productType = esId,
                    .reactionRadius = 0.03,
                    .rate = 86.78638438
            };

        }
    }


    using ExternalPotentials = std::tuple<>;
    using PairPotentials = std::tuple<>;

    using ReactionsO1 = std::tuple<ctiprd::reactions::doi::Fission<T>, ctiprd::reactions::doi::Fission<T>>;
    using ReactionsO2 = std::tuple<ctiprd::reactions::doi::Fusion<T>>;

    ReactionsO1 reactionsO1{};
    ReactionsO2 reactionsO2{};

    ExternalPotentials externalPotentials{};
    PairPotentials pairPotentials{};
};
}
