//
// Created by mho on 3/28/22.
//

#include <catch2/catch.hpp>
#include <ctiprd/cpu/UncontrolledApproximation.h>
#include <ctiprd/systems/double_well.h>

TEST_CASE("UncontrolledApproximation sanity", "[reactions]") {
    using System = ctiprd::systems::DoubleWell<float>;
    using ParticleCollection = ctiprd::cpu::ParticleCollection<System, ctiprd::cpu::particles::positions>;
    ctiprd::cpu::UncontrolledApproximation<ParticleCollection, System> ua {System{}};
}
