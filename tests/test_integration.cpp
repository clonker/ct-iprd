#include <catch2/catch.hpp>

#include <ctiprd/NeighborList.h>
#include "ctiprd/systems/double_well.h"
#include "spdlog/spdlog.h"

using System = ctiprd::systems::DoubleWell<float>;

TEST_CASE("Test Integration", "[integration]") {
    std::size_t nSteps = 100;

    System system {};
    auto pool = ctiprd::config::make_pool(3);
    auto integrator = ctiprd::integrator::EulerMaruyama{system, pool};
    for(int n = 0; n < 500; ++n) {
        integrator.particles()->addParticle({{0., 0.}}, "A");
    }
    for(int n = 0; n < 500; ++n) {
        integrator.particles()->addParticle({{.1, .1}}, "B");
    }


    spdlog::set_level(spdlog::level::debug);
    for(std::size_t t = 0; t < nSteps; ++t) {
        integrator.step(1e-2);

        spdlog::debug("{}: n particles {}", t, integrator.particles()->nParticles());
    }
}
