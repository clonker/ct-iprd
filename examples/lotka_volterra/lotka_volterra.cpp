//
// Created by mho on 4/5/22.
//

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

    m.def("simulate", [] (std::size_t nSteps, float dt, int njobs, const std::function<void(std::size_t)>& progressCallback) {
        System system {};
        auto pool = ctiprd::config::make_pool(njobs);
        auto integrator = ctiprd::integrator::EulerMaruyama{system, pool};

        {
            auto &generator = ctiprd::rnd::staticThreadLocalGenerator();
            std::uniform_real_distribution<float> d1 {-System::boxSize[0] / 2, System::boxSize[0] / 2};
            std::uniform_real_distribution<float> d2 {-System::boxSize[1] / 2, System::boxSize[1] / 2};
            for (int n = 0; n < 125; ++n) {
                integrator.particles()->addParticle({{d1(generator), d2(generator)}}, "prey");
            }
            for (int n = 0; n < 100; ++n) {
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
                    trajectory.emplace_back(
                            np_array<float, nb::any, 2>{},
                            np_array<std::size_t, nb::any>{});
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
