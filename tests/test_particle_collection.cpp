#include <set>

#include <catch2/catch.hpp>

#include <ctiprd/systems/double_well.h>
#include <ctiprd/systems/michaelis_menten.h>
#include <ctiprd/cpu/integrators/EulerMaruyama.h>

TEST_CASE("Particle collection sanity", "[particles]") {
    using System = ctiprd::systems::DoubleWell<float>;
    ctiprd::cpu::ParticleCollection<System, ctiprd::cpu::particles::positions, ctiprd::cpu::particles::forces> collection {};
    REQUIRE(collection.dim == 2);
    REQUIRE(collection.containsForces() == true);
    REQUIRE(collection.containsPositions() == true);
    REQUIRE(collection.containsVelocities() == false);
}

SCENARIO("Particle collection can have different flags") {
    using System = ctiprd::systems::DoubleWell<float>;
    auto pool = ctiprd::config::make_pool(5);
    using Collection = ctiprd::cpu::ParticleCollection<System, ctiprd::cpu::particles::positions, ctiprd::cpu::particles::forces>;
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
                pool->waitForTasks();

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
        ctiprd::cpu::ParticleCollection<System, ctiprd::cpu::particles::positions> collection {};
        THEN("the flags are set up correctly") {
            REQUIRE(collection.containsForces() == false);
            REQUIRE(collection.containsPositions() == true);
            REQUIRE(collection.containsVelocities() == false);
        }
    }

    GIVEN("A particle collection with positions and velocities") {
        ctiprd::cpu::ParticleCollection<System, ctiprd::cpu::particles::positions, ctiprd::cpu::particles::velocities> collection {};
        THEN("the flags are set up correctly") {
            REQUIRE(collection.containsForces() == false);
            REQUIRE(collection.containsPositions() == true);
            REQUIRE(collection.containsVelocities() == true);
        }
    }
}

TEST_CASE("Test against michaelis menten", "[particles][michaelis-menten]") {
    using System = ctiprd::systems::MichaelisMenten<float>;
    System system {};
    auto pool = ctiprd::config::make_pool(8);
    auto integrator = ctiprd::cpu::integrator::EulerMaruyama{system, pool};
    int nEinit = 90;
    int nSinit = 90;
    {
        auto &generator = ctiprd::rnd::staticThreadLocalGenerator();
        std::uniform_real_distribution<float> d1 {-System::boxSize[0] / 2, System::boxSize[0] / 2};
        std::uniform_real_distribution<float> d2 {-System::boxSize[1] / 2, System::boxSize[1] / 2};
        std::uniform_real_distribution<float> d3 {-System::boxSize[2] / 2, System::boxSize[2] / 2};
        for (int n = 0; n < nEinit; ++n) {
            integrator.particles()->addParticle({{d1(generator), d2(generator), d3(generator)}}, "E");
        }
        for (int n = 0; n < nSinit; ++n) {
            integrator.particles()->addParticle({{d1(generator), d2(generator), d3(generator)}}, "S");
        }
    }

    REQUIRE(integrator.particles()->nParticles() == nEinit + nSinit);
    std::atomic<std::size_t> nE {0};
    std::atomic<std::size_t> nS {0};
    std::atomic<std::size_t> nES {0};
    std::atomic<std::size_t> nP {0};
    {
        integrator.particles()->forEachParticle([&nE, &nS, &nES, &nP](const auto &id, const auto &, const auto &type,
                                                                          const auto &) {

            if(type == System::eId) {
                ++nE;
            } else if(type == System::sId) {
                ++nS;
            } else if(type == System::esId) {
                ++nES;
            } else if(type == System::pId) {
                ++nP;
            }
        }, pool);
        pool->waitForTasks();
    }
    REQUIRE(nE.load() == nEinit);
    REQUIRE(nS.load() == nSinit);
    REQUIRE(nES.load() == 0);
    REQUIRE(nP.load() == 0);
}
