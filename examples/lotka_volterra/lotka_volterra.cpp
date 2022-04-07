//
// Created by mho on 4/5/22.
//

#include <memory>

#include <ctiprd/systems/lotka_volterra.h>
#include <ctiprd/binding/system_bindings.h>
#include <ctiprd/progressbar.hpp>

namespace nb = nanobind;

template<typename T, std::size_t... shape>
using np_array = ctiprd::binding::np_array<T, shape...>;

using System = ctiprd::systems::LotkaVolterra<float>;

NB_MODULE(lv_mod, m) {
    ctiprd::binding::exportBaseTypes<System::dtype>(m);
    ctiprd::binding::exportSystem<System>(m, "LotkaVolterra");

    m.def("simulate", [] (std::size_t nSteps, float dt, int njobs, nb::handle progressCallback) {
        System system {};
        auto pool = ctiprd::config::make_pool(njobs);
        auto integrator = ctiprd::integrator::EulerMaruyama{system, pool};

        {
            auto &generator = ctiprd::rnd::staticThreadLocalGenerator();
            std::uniform_real_distribution<float> d1 {-System::boxSize[0] / 2, System::boxSize[0] / 2};
            std::uniform_real_distribution<float> d2 {-System::boxSize[1] / 2, System::boxSize[1] / 2};
            for (int n = 0; n < 9600; ++n) {
                integrator.particles()->addParticle({{d1(generator), d2(generator)}}, "prey");
            }
            for (int n = 0; n < 9600; ++n) {
                integrator.particles()->addParticle({{d1(generator), d2(generator)}}, "predator");
            }
        }

        std::vector<std::tuple<np_array<float, nb::any, 2>, np_array<std::size_t, nb::any>>> trajectory;

        {
            nb::gil_scoped_release release;
            for (std::size_t t = 0; t < nSteps; ++t) {
                integrator.step(dt);

                if (t % 200 == 0) {
                    nb::gil_scoped_acquire acquire;

                    auto nParticles = integrator.particles()->nParticles();
                    auto* ptrTraj = new float[2 * nParticles];
                    auto* ptrTypes = new std::size_t[nParticles];
                    nb::capsule deleterTraj(ptrTraj, [](void *data) {
                        delete[] (float *) data;
                    });
                    nb::capsule deleterTypes(ptrTypes, [](void *data) {
                        delete[] (std::size_t *) data;
                    });

                    std::size_t shapeTraj[2] = {nParticles, 2};
                    std::size_t shapeTypes[1] = {nParticles};
                    trajectory.emplace_back(
                            np_array<float, nb::any, 2>{ptrTraj, 2, shapeTraj, deleterTraj},
                            np_array<std::size_t, nb::any>{ptrTypes, 1, shapeTypes, deleterTypes});
                    std::size_t nPredator{}, nPrey{};
                    std::size_t ix = 0;
                    for (std::size_t i = 0; i < integrator.particles()->size(); ++i) {
                        if (integrator.particles()->exists(i)) {
                            std::get<0>(trajectory.back())(ix, 0) = integrator.particles()->positionOf(i)[0];
                            std::get<0>(trajectory.back())(ix, 1) = integrator.particles()->positionOf(i)[1];
                            std::get<1>(trajectory.back())(ix) = integrator.particles()->typeOf(i);
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
                    nb::gil_scoped_acquire acquire;
                    if (PyErr_CheckSignals() != 0) {
                        throw nb::python_error();
                    }
                }
            }
        }

        return trajectory;
    });
}
