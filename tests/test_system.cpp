//
// Created by mho on 3/25/22.
//

#include <catch2/catch.hpp>
#include <ctiprd/systems/double_well.h>

using System = ctiprd::systems::DoubleWell<float>;

TEST_CASE("System sanity", "[system]") {
    REQUIRE(std::tuple_element_t<0, System::ExternalPotentials>::particleType == 0);
    REQUIRE(std::tuple_element_t<1, System::ExternalPotentials>::particleType == 1);
    REQUIRE(std::tuple_element_t<0, System::PairPotentials>::particleType1 == 0);
    REQUIRE(std::tuple_element_t<0, System::PairPotentials>::particleType2 == 0);
    REQUIRE(std::tuple_element_t<1, System::PairPotentials>::particleType1 == 0);
    REQUIRE(std::tuple_element_t<1, System::PairPotentials>::particleType2 == 1);
    REQUIRE(std::tuple_element_t<2, System::PairPotentials>::particleType1 == 1);
    REQUIRE(std::tuple_element_t<2, System::PairPotentials>::particleType2 == 1);
}
