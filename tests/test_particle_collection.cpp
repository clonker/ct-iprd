#include <catch2/catch.hpp>
#include <ctiprd/systems/double_well.h>
#include "ctiprd/ParticleCollection.h"

using System = ctiprd::systems::DoubleWell<float>;

TEST_CASE("Particle collection sanity", "[particles]") {
    ctiprd::ParticleCollection<System> collection {};
    REQUIRE(collection.dim == 2);
    REQUIRE(collection.containsForces() == true);
    REQUIRE(collection.containsPositions() == true);
    REQUIRE(collection.containsVelocities() == false);
}

SCENARIO("Particle collection can have different flags") {
    auto pool = ctiprd::config::make_pool(5);
    using Collection = ctiprd::ParticleCollection<System>;
    GIVEN("A particle collection with default flags") {
        Collection collection {};
        THEN("the flags are set up correctly") {
            REQUIRE(collection.containsForces() == true);
            REQUIRE(collection.containsPositions() == true);
            REQUIRE(collection.containsVelocities() == false);
        }

        WHEN("Adding a particle") {
            collection.addParticle({{0., 0.}}, "A");
            THEN("The size grows") {
                REQUIRE(collection.nParticles() == 1);
            }

        }

        WHEN("Adding multiple particles") {
            for(int i = 0; i < 1000; ++i) collection.addParticle({{0., 0.}}, "A");

            THEN("The size grows") {
                REQUIRE(collection.nParticles() == 1000);
            }

            WHEN("Modifying the particles") {
                auto futures = collection.forEachParticle([](std::size_t i, Collection::Position &p, auto type, Collection::Force &f) {
                    p[0] = 55;
                    p[1] = 22;

                    f[0] = 11;
                    f[1] = -11;
                }, pool);
                std::for_each(begin(futures), end(futures), [](auto &f) { f.wait(); });

                THEN("This modification is reflected in the data") {
                    Collection::ContainerType<Collection::MaybePosition> vecRef (1000, {{55, 22}});
                    REQUIRE(collection.positions() == vecRef);

                    Collection::ContainerType<Collection::Force> vecRefForces (1000, {{11, -11}});
                    REQUIRE(collection.forces() == vecRefForces);
                }
            }
        }
    }

    GIVEN("A particle collection with just positions") {
        ctiprd::ParticleCollection<System, ctiprd::particle_collection::usePositions> collection {};
        THEN("the flags are set up correctly") {
            REQUIRE(collection.containsForces() == false);
            REQUIRE(collection.containsPositions() == true);
            REQUIRE(collection.containsVelocities() == false);
        }
    }

    GIVEN("A particle collection with positions and velocities") {
        ctiprd::ParticleCollection<System, ctiprd::particle_collection::usePositions | ctiprd::particle_collection::useVelocities> collection {};
        THEN("the flags are set up correctly") {
            REQUIRE(collection.containsForces() == false);
            REQUIRE(collection.containsPositions() == true);
            REQUIRE(collection.containsVelocities() == true);
        }
    }
}
