#include <cstddef>
#include <tuple>
#include <optional>

#include <fmt/format.h>

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <ctiprd/systems/michaelis_menten.h>

#include <ctiprd/reactions/doi.h>
#include <ctiprd/systems/util.h>
#include <ctiprd/binding/system_bindings.h>
#include <ctiprd/cpu/integrators/EulerMaruyama.h>

namespace py = pybind11;

using System = ctiprd::systems::MichaelisMenten<float>;

PYBIND11_MODULE(mm_mod, m) {
    ctiprd::binding::exportBaseTypes<System::dtype>(m);
    ctiprd::binding::exportSystem<System>(m, "MichaelisMenten");

    m.def("simulate", [] (std::size_t nSteps, float dt, int njobs, py::handle progressCallback) {
        System system {};
        auto pool = ctiprd::config::make_pool(njobs);
        auto integrator = ctiprd::cpu::integrator::EulerMaruyama{system, pool};

        integrator.particles()->initializeParticles(909, "E");
        integrator.particles()->initializeParticles(9091, "S");

        std::vector<std::size_t> nE {};
        std::vector<std::size_t> nS {};
        std::vector<std::size_t> nES {};
        std::vector<std::size_t> nP {};
        {
            py::gil_scoped_release release;
            for (std::size_t t = 0; t < nSteps; ++t) {
                {
                    std::atomic<std::size_t> nEl {0};
                    std::atomic<std::size_t> nSl {0};
                    std::atomic<std::size_t> nESl {0};
                    std::atomic<std::size_t> nPl {0};
                    auto futures = integrator.particles()->forEachParticle([&nEl, &nSl, &nESl, &nPl](const auto &id, const auto &, const auto &type,
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
                    for(auto &f : futures) f.wait();

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
