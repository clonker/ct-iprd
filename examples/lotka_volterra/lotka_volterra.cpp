//
// Created by mho on 4/5/22.
//

#include <memory>

#include <ctiprd/systems/lotka_volterra.h>
#include <ctiprd/binding/system_bindings.h>
#include <ctiprd/progressbar.hpp>

namespace py = pybind11;

template<typename T>
using np_array = ctiprd::binding::np_array<T>;
using System = ctiprd::systems::LotkaVolterra<float>;

PYBIND11_MODULE(lv_mod, m) {
    ctiprd::binding::exportBaseTypes<System::dtype>(m);
    ctiprd::binding::exportSystem<System>(m, "LotkaVolterra");

    m.def("simulate", [] (std::size_t nSteps, float dt, int njobs, int nPrey, int nPredator, py::handle progressCallback) {
        System system {};
        auto pool = ctiprd::config::make_pool(njobs);
        auto integrator = ctiprd::integrator::EulerMaruyama{system, pool};

        {
            auto &generator = ctiprd::rnd::staticThreadLocalGenerator();
            std::uniform_real_distribution<float> d1 {-System::boxSize[0] / 2, System::boxSize[0] / 2};
            std::uniform_real_distribution<float> d2 {-System::boxSize[1] / 2, System::boxSize[1] / 2};
            for (int n = 0; n < nPrey; ++n) {
                integrator.particles()->addParticle({{d1(generator), d2(generator)}}, "prey");
            }
            for (int n = 0; n < nPredator; ++n) {
                integrator.particles()->addParticle({{d1(generator), d2(generator)}}, "predator");
            }
        }

        std::vector<std::tuple<np_array<float>, np_array<std::size_t>>> trajectory;

        {
            py::gil_scoped_release release;
            for (std::size_t t = 0; t < nSteps; ++t) {
                integrator.step(dt);

                if (t % 200 == 0) {
                    py::gil_scoped_acquire acquire;

                    auto nParticles = integrator.particles()->nParticles();
                    auto* ptrTraj = new float[2 * nParticles];
                    auto* ptrTypes = new std::size_t[nParticles];
                    py::capsule deleterTraj(ptrTraj, [](void *data) {
                        delete[] (float *) data;
                    });
                    py::capsule deleterTypes(ptrTypes, [](void *data) {
                        delete[] (std::size_t *) data;
                    });

                    std::size_t shapeTraj[2] = {nParticles, 2};
                    std::size_t shapeTypes[1] = {nParticles};
                    trajectory.emplace_back(np_array<float>{shapeTraj}, np_array<std::size_t>{shapeTypes});
                    std::size_t ix = 0;
                    for (std::size_t i = 0; i < integrator.particles()->size(); ++i) {
                        if (integrator.particles()->exists(i)) {
                            std::get<0>(trajectory.back()).mutable_at(ix, 0) = integrator.particles()->positionOf(i)[0];
                            std::get<0>(trajectory.back()).mutable_at(ix, 1) = integrator.particles()->positionOf(i)[1];
                            std::get<1>(trajectory.back()).mutable_at(ix) = integrator.particles()->typeOf(i);
                            if (integrator.particles()->typeOf(i) == System::preyId) {
                                ++nPrey;
                            } else {
                                ++nPredator;
                            }
                            ++ix;
                        }
                    }
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

        return trajectory;
    });
}
