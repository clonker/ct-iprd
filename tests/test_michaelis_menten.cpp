#include <catch2/catch.hpp>

#include <ctiprd/systems/michaelis_menten.h>

#include <ctiprd/reactions/doi.h>
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

            {
                std::atomic<std::size_t> nEl {0};
                std::atomic<std::size_t> nSl {0};
                std::atomic<std::size_t> nESl {0};
                std::atomic<std::size_t> nPl {0};
                integrator.particles()->forEachParticle([&nEl, &nSl, &nESl, &nPl](const auto &id, const auto &, const auto &type,
                                                                                  const auto &) {
                    if(type == System::eId) {
                        ++nEl;
                    } else if(type == System::sId) {
                        ++nSl;
                    } else if(type == System::esId) {
                        ++nESl;
                    } else if(type == System::pId) {
                        ++nPl;
                    }
                }, pool);
                pool->waitForTasks();

                nE.emplace_back(nEl.load());
                nS.emplace_back(nSl.load());
                nES.emplace_back(nESl.load());
                nP.emplace_back(nPl.load());
            }

            integrator.step(8e-4);
        }
    }
}

