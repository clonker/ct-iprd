#include <catch2/catch.hpp>
#include "ctiprd/ParticleCollection.h"

TEST_CASE("Particle collection sanity", "[particles]") {
    ctiprd::ParticleCollection<3, float> collection {"A"};
    REQUIRE(collection.dim == 3);
    REQUIRE(collection.containsForces() == true);
    REQUIRE(collection.containsPositions() == true);
    REQUIRE(collection.containsVelocities() == true);
}

SCENARIO("Particle collection can have different flags") {
    GIVEN("A particle collection with all flags") {
        ctiprd::ParticleCollection<3, float> collection {"A"};
        THEN("the flags are set up correctly") {
            REQUIRE(collection.containsForces() == true);
            REQUIRE(collection.containsPositions() == true);
            REQUIRE(collection.containsVelocities() == true);
        }

        WHEN("Adding a particle") {
            collection.addParticle({0, 0, 0}, {0, 0, 0}, {0, 0, 0});
        }
    }

    GIVEN("A particle collection with just positions") {
        ctiprd::ParticleCollection<3, float, ctiprd::particle_collection::Flags::usePositions> collection {"A"};
        THEN("the flags are set up correctly") {
            REQUIRE(collection.containsForces() == false);
            REQUIRE(collection.containsPositions() == true);
            REQUIRE(collection.containsVelocities() == false);
        }
    }

    GIVEN("A particle collection with positions and velocities") {
        ctiprd::ParticleCollection<3, float, ctiprd::particle_collection::Flags::usePositions | ctiprd::particle_collection::Flags::useVelocities> collection {"A"};
        THEN("the flags are set up correctly") {
            REQUIRE(collection.containsForces() == false);
            REQUIRE(collection.containsPositions() == true);
            REQUIRE(collection.containsVelocities() == true);
        }
    }
}
