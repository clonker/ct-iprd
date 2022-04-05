//
// Created by mho on 4/5/22.
//

#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <ctiprd/systems/lotka_volterra.h>
#include <ctiprd/progressbar.hpp>

namespace py = pybind11;

using System = ctiprd::systems::LotkaVolterra<float>;

template<typename dtype>
using np_array = py::array_t<dtype, py::array::c_style | py::array::forcecast>;

PYBIND11_MODULE(lv_mod, m) {
    m.def("simulate", [] (std::size_t nSteps, int njobs, std::function<void(std::size_t)> progressCallback) {
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

        std::vector<std::tuple<np_array<float>, np_array<std::size_t>>> trajectory;

        // progressbar bar(static_cast<int>(nSteps));
        for(std::size_t t = 0; t < nSteps; ++t) {
            py::gil_scoped_release gilScopedRelease;
            integrator.step(1e-2);

            if(t % 200 == 0) {
                py::gil_scoped_acquire gil;
                trajectory.emplace_back(
                        np_array<float>{std::vector<std::size_t>{integrator.particles()->nParticles(), 2}},
                        np_array<std::size_t>{std::vector<std::size_t>(1, integrator.particles()->nParticles())});
                std::size_t nPredator {}, nPrey {};
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
                // spdlog::critical("nPrey = {}, nPredator = {}", nPrey, nPredator);
            }
            // bar.update();

            if (PyErr_CheckSignals() != 0){
                throw py::error_already_set();
            }
        }

        return trajectory;
    });
}
