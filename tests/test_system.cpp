//
// Created by mho on 3/25/22.
//

#include <catch2/catch.hpp>
#include <ctiprd/systems/double_well.h>

using System = ctiprd::systems::DoubleWell<float>;

TEST_CASE("System sanity", "[system]") {
    System system;
    auto &[pot1, pot2] = system.externalPotentials;
    REQUIRE(pot1.particleType == 0);
    REQUIRE(pot2.particleType == 1);
}
