#include <catch2/catch.hpp>

#include <spdlog/spdlog.h>

#include <ctiprd/systems/double_well.h>
#include <ctiprd/systems/lotka_volterra.h>

#include <ctiprd/cpu/integrators/EulerMaruyama.h>

TEMPLATE_TEST_CASE("Test Integration", "[integration]", ctiprd::systems::DoubleWell<float>, ctiprd::systems::LotkaVolterra<float>) {
    std::size_t nSteps = 100;

    TestType system {};
    auto pool = ctiprd::config::make_pool(6);
    auto integrator = ctiprd::cpu::integrator::EulerMaruyama{system, pool};

    {
        auto &generator = ctiprd::rnd::staticThreadLocalGenerator();
        for(std::size_t i = 0; i < 10000; ++i) {
            ctiprd::Vec<float, 2> pos {};
            for (std::size_t d = 0; d < 2; ++d) {
                std::uniform_real_distribution<float> dist {-TestType::boxSize[d], TestType::boxSize[d]};
                pos[d] = dist(generator);
            }

            integrator.particles()->addParticle(pos, 0);
        }
    }
    integrator.particles()->initializeParticles(25000, system.types[1].name);
    /*integrator.particles()->initializeParticles(25000, system.types[0].name);
    integrator.particles()->initializeParticles(25000, system.types[1].name);*/

    spdlog::set_level(spdlog::level::debug);
    for(std::size_t t = 0; t < nSteps; ++t) {
        integrator.step(1e-4);

        spdlog::debug("{}: n particles {}", t, integrator.particles()->nParticles());
    }
}
