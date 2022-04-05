//
// Created by mho on 4/5/22.
//

#include <pybind11/numpy.h>

#include <ctiprd/systems/lotka_volterra.h>
#include <ctiprd/progressbar.hpp>

namespace py = pybind11;

using System = ctiprd::systems::LotkaVolterra<float>;

template<typename dtype>
using np_array = py::array_t<dtype, py::array::c_style | py::array::forcecast>;

PYBIND11_MODULE(lv_mod, m) {
    m.def("simulate", [] (int njobs) {
        std::size_t nSteps = 100000;

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

        // np_array<float> out {{nSteps, static_cast<std::size_t>(integrator.particles()->nParticles()), static_cast<std::size_t>(2)}};

        progressbar bar(static_cast<int>(nSteps));
        for(std::size_t t = 0; t < nSteps; ++t) {
            integrator.step(1e-2);

            /*const auto &positions = integrator.particles()->positions();
            auto i = 0;
            for (const auto &pos: positions) {
                out.mutable_at(t, i, 0) = (*pos)[0];
                out.mutable_at(t, i, 1) = (*pos)[1];
                ++i;
            }*/

            bar.update();
        }

        // return out;
    });
}
