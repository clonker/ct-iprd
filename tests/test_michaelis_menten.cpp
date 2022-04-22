#include <catch2/catch.hpp>

#include <ctiprd/systems/michaelis_menten.h>

#include <ctiprd/systems/util.h>
#include <ctiprd/cpu/integrators/EulerMaruyama.h>


TEST_CASE("Michaelis Menten", "[michaelis-menten][integration]") {
    using System = ctiprd::systems::MichaelisMenten<float>;

    System system {};
    auto pool = ctiprd::config::make_pool(5);
    auto integrator = ctiprd::cpu::integrator::EulerMaruyama{system, pool};

    integrator.particles()->initializeParticles(909, "E");
    integrator.particles()->initializeParticles(9091, "S");

    std::vector<std::size_t> nE {};
    std::vector<std::size_t> nS {};
    std::vector<std::size_t> nES {};
    std::vector<std::size_t> nP {};

    std::size_t nSteps {100};
    {
        for (std::size_t t = 0; t < nSteps; ++t) {
            integrator.step(8e-4);
        }
    }
}

