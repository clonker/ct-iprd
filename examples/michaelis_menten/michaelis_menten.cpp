#include <cstddef>
#include <tuple>
#include <optional>

#include <fmt/format.h>

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <ctiprd/systems/michaelis_menten.h>

#include <ctiprd/integrators/EulerMaruyama.h>
#include <ctiprd/reactions/basic.h>
#include <ctiprd/systems/util.h>
#include <ctiprd/binding/system_bindings.h>

namespace py = pybind11;



using System = ctiprd::systems::MichaelisMenten<float>;

PYBIND11_MODULE(mm_mod, m) {
    ctiprd::binding::exportBaseTypes<System::dtype>(m);
    ctiprd::binding::exportSystem<System>(m, "MichaelisMenten");

    m.def("simulate", [] (std::size_t nSteps, float dt, int njobs, py::handle progressCallback) {
        System system {};
        auto pool = ctiprd::config::make_pool(njobs);
        auto integrator = ctiprd::integrator::EulerMaruyama{system, pool};

        {
            auto &generator = ctiprd::rnd::staticThreadLocalGenerator();
            std::uniform_real_distribution<float> d1 {-System::boxSize[0] / 2, System::boxSize[0] / 2};
            std::uniform_real_distribution<float> d2 {-System::boxSize[1] / 2, System::boxSize[1] / 2};
            std::uniform_real_distribution<float> d3 {-System::boxSize[2] / 2, System::boxSize[2] / 2};
            int nE = 909;
            int nS = 9091;
            for (int n = 0; n < nE; ++n) {
                integrator.particles()->addParticle({{d1(generator), d2(generator), d3(generator)}}, "E");
            }
            for (int n = 0; n < nS; ++n) {
                integrator.particles()->addParticle({{d1(generator), d2(generator), d3(generator)}}, "S");
            }
        }
        std::vector<std::size_t> nE {};
        std::vector<std::size_t> nS {};
        std::vector<std::size_t> nES {};
        std::vector<std::size_t> nP {};
        std::mutex m {};
        {
            py::gil_scoped_release release;
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

                integrator.step(dt);
                if (t % 200 == 0) {
                    py::gil_scoped_acquire acquire;
                    progressCallback(t);
                }

                if(t % 50 == 0) {
                    py::gil_scoped_acquire acquire;
                    if (PyErr_CheckSignals() != 0) {
                        throw py::error_already_set();
                    }
                }
            }
        }

        return std::make_tuple(nE, nS, nES, nP);
    });
}
