//
// Created by mho on 4/13/22.
//

#include <catch2/catch.hpp>
#include <ctiprd/systems/util.h>
#include <ctiprd/cpu/reactions.h>
#include "ctiprd/ParticleTypes.h"
#include "ctiprd/cpu/ParticleCollection.h"

template<typename T>
struct TestSystem {
    using dtype = T;
    static constexpr std::size_t DIM = 2;
    static constexpr std::array<T, DIM> boxSize {5., 5.};
    static constexpr bool periodic = true;
    static constexpr T kBT = 1.;
    static constexpr ctiprd::ParticleTypes<dtype, 3> types {{
                                                                    {
                                                                            .name = "A",
                                                                            .diffusionConstant = 1.
                                                                    },
                                                                    {
                                                                            .name = "B",
                                                                            .diffusionConstant = 1.
                                                                    },
                                                            }};
    static constexpr auto aId = ctiprd::systems::particleTypeId<types>("A");
    static constexpr auto bId = ctiprd::systems::particleTypeId<types>("B");

    TestSystem() {
        {
            auto &[conv1, conv2, decay, fission] = reactionsO1;
            conv1.eductType = aId;
            conv1.productType = bId;
            conv1.rate = 1.;

            conv2.eductType = bId;
            conv2.productType = aId;
            conv2.rate = 2.;

            decay.eductType = aId;
            decay.rate = 3.;

            fission.eductType = bId;
            fission.rate = 4.;
            fission.productType1 = aId;
            fission.productType2 = aId;
            fission.productDistance = 2.;
        }
        {
            auto &[fusion, cata] = reactionsO2;
            fusion.eductType1 = aId;
            fusion.eductType2 = aId;
            fusion.productType = bId;
            fusion.rate = .1;
            fusion.reactionRadius = .1;

            cata.eductType = aId;
            cata.catalyst = bId;
            cata.productType = bId;
            cata.rate = .1;
            cata.reactionRadius = .1;
        }
    }



    using ExternalPotentials = std::tuple<>;
    using PairPotentials = std::tuple<>;
    using ReactionsO1 = std::tuple<ctiprd::reactions::doi::Conversion<T>,
                                   ctiprd::reactions::doi::Conversion<T>,
                                   ctiprd::reactions::doi::Decay<T>,
                                   ctiprd::reactions::doi::Fission<T>>;
    using ReactionsO2 = std::tuple<ctiprd::reactions::doi::Fusion<T>, ctiprd::reactions::doi::Catalysis<T>>;

    ReactionsO1 reactionsO1 {};
    ReactionsO2 reactionsO2 {};

    ExternalPotentials externalPotentials {};
    PairPotentials pairPotentials {};
};


TEST_CASE("Test reactions map", "[reactions]") {
    using System = TestSystem<float>;
    using ParticleCollection = ctiprd::cpu::ParticleCollection<System, ctiprd::cpu::particles::positions>;
    using Updater = ctiprd::cpu::ParticleCollectionUpdater<System, ParticleCollection>;
    System system {};
    {
        auto[map, backing] = ctiprd::cpu::reactions::impl::generateMapO1<Updater>(system);
        REQUIRE(map[System::aId].size() == 2);
        REQUIRE(map[System::aId][0]->shouldPerform(1e50, {}, System::aId));
        REQUIRE(!map[System::aId][0]->shouldPerform(1e50, {}, System::bId));
        REQUIRE(map[System::aId][1]->shouldPerform(1e50, {}, System::aId));
        REQUIRE(!map[System::aId][1]->shouldPerform(1e50, {}, System::bId));

        REQUIRE(map[System::bId].size() == 2);
        REQUIRE(map[System::bId][0]->shouldPerform(1e50, {}, System::bId));
        REQUIRE(!map[System::bId][0]->shouldPerform(1e50, {}, System::aId));
        REQUIRE(map[System::bId][1]->shouldPerform(1e50, {}, System::bId));
        REQUIRE(!map[System::bId][1]->shouldPerform(1e50, {}, System::aId));
    }
    {
        auto [map, backing] = ctiprd::cpu::reactions::impl::generateMapO2<Updater>(system);
        REQUIRE(map[std::make_tuple(System::aId, System::aId)].size() == 1);
        REQUIRE(!map[std::make_tuple(System::aId, System::aId)][0]->shouldPerform(1e50, {}, System::aId, {}, System::bId));
        REQUIRE(!map[std::make_tuple(System::aId, System::aId)][0]->shouldPerform(1e50, {}, System::bId, {}, System::aId));
        REQUIRE(!map[std::make_tuple(System::aId, System::aId)][0]->shouldPerform(1e50, {}, System::bId, {}, System::bId));
        REQUIRE(map[std::make_tuple(System::aId, System::aId)][0]->shouldPerform(1e50, {}, System::aId, {}, System::aId));

        REQUIRE(map[std::make_tuple(System::aId, System::bId)].size() == 1);
        REQUIRE(map[std::make_tuple(System::bId, System::aId)].size() == 1);
        REQUIRE(map[std::make_tuple(System::bId, System::aId)] == map[std::make_tuple(System::aId, System::bId)]);
        REQUIRE(map[std::make_tuple(System::aId, System::bId)][0]->shouldPerform(1e50, {}, System::aId, {}, System::bId));
        REQUIRE(map[std::make_tuple(System::aId, System::bId)][0]->shouldPerform(1e50, {}, System::bId, {}, System::aId));
        REQUIRE(!map[std::make_tuple(System::aId, System::bId)][0]->shouldPerform(1e50, {}, System::aId, {}, System::aId));
        REQUIRE(!map[std::make_tuple(System::aId, System::bId)][0]->shouldPerform(1e50, {}, System::bId, {}, System::bId));
    }
}
