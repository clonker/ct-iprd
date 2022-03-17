#include <catch2/catch.hpp>
#include "ctiprd/ParticleCollection.h"

TEST_CASE("Particle collection sanity", "[particles]") {
    ctiprd::ParticleCollection<3, float> collection {"A"};
    REQUIRE(collection.dim == 3);
    REQUIRE(collection.containsForces() == true);
    REQUIRE(collection.containsPositions() == true);
    REQUIRE(collection.containsVelocities() == false);
}

SCENARIO("Particle collection can have different flags") {
    auto pool = ctiprd::config::make_pool(5);
    using Collection = ctiprd::ParticleCollection<2, float>;
    GIVEN("A particle collection with default flags") {
        Collection collection {"A"};
        THEN("the flags are set up correctly") {
            REQUIRE(collection.containsForces() == true);
            REQUIRE(collection.containsPositions() == true);
            REQUIRE(collection.containsVelocities() == false);
        }

        WHEN("Adding a particle") {
            collection.addParticle({{0., 0.}});
            THEN("The size grows") {
                REQUIRE(collection.nParticles() == 1);
            }

        }

        WHEN("Adding multiple particles") {
            for(int i = 0; i < 1000; ++i) collection.addParticle({{0., 0.}});

            THEN("The size grows") {
                REQUIRE(collection.nParticles() == 1000);
            }

            WHEN("Modifying the particles") {
                auto futures = collection.forEachParticle([](std::size_t i, Collection::Position &p, Collection::Force &f) {
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
        ctiprd::ParticleCollection<3, float, ctiprd::particle_collection::usePositions> collection {"A"};
        THEN("the flags are set up correctly") {
            REQUIRE(collection.containsForces() == false);
            REQUIRE(collection.containsPositions() == true);
            REQUIRE(collection.containsVelocities() == false);
        }
    }

    GIVEN("A particle collection with positions and velocities") {
        ctiprd::ParticleCollection<3, float, ctiprd::particle_collection::usePositions | ctiprd::particle_collection::useVelocities> collection {"A"};
        THEN("the flags are set up correctly") {
            REQUIRE(collection.containsForces() == false);
            REQUIRE(collection.containsPositions() == true);
            REQUIRE(collection.containsVelocities() == true);
        }
    }
}
