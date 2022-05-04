//
// Created by mho on 3/16/22.
//

#include <pybind11/numpy.h>

#include <ctiprd/config.h>

#include <ctiprd/systems/double_well.h>
#include <ctiprd/progressbar.hpp>
#include <ctiprd/cpu/integrators/EulerMaruyama.h>

namespace py = pybind11;

using System = ctiprd::systems::DoubleWell<float>;

template<typename dtype>
using np_array = py::array_t<dtype, py::array::c_style | py::array::forcecast>;

PYBIND11_MODULE(dw_mod, m) {
    m.def("simulate", [] () {
        std::size_t nSteps = 1000;

        System system {};
        auto pool = ctiprd::config::make_pool(3);
        auto integrator = ctiprd::cpu::integrator::EulerMaruyama{system, pool};
        for(int n = 0; n < 10000; ++n) {
            integrator.particles()->addParticle({{0., 0.}}, "A");
        }

        np_array<float> out {{nSteps, static_cast<std::size_t>(integrator.particles()->nParticles()), static_cast<std::size_t>(2)}};

        progressbar bar(nSteps);
        for(std::size_t t = 0; t < nSteps; ++t) {
            integrator.step(1e-3);

            const auto &positions = integrator.particles()->positions_();
            auto i = 0;
            for (const auto &pos: positions) {
                out.mutable_at(t, i, 0) = (*pos)[0];
                out.mutable_at(t, i, 1) = (*pos)[1];
                ++i;
            }

            bar.update();
        }

        return out;
    });
}
