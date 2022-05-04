//
// Created by mho on 4/5/22.
//

#include <memory>

#include <pybind11/stl.h>

#include <ctiprd/config.h>

#include <ctiprd/systems/lotka_volterra.h>
#include <ctiprd/binding/system_bindings.h>
#include <ctiprd/progressbar.hpp>

#include <ctiprd/cpu/integrators/EulerMaruyama.h>

namespace py = pybind11;

template<typename T>
using np_array = ctiprd::binding::np_array<T>;
using System = ctiprd::systems::LotkaVolterra<double>;

void check_shape(const np_array<System::dtype> &arr) {
    if(arr.ndim() != 2) {
        throw std::runtime_error(
                fmt::format("Provided particle array needs to be two-dimensional but was {}-dimensional.", arr.ndim())
        );
    }
    if(arr.shape(1) != System::DIM) {
        throw std::runtime_error(
                fmt::format("System is {}-dimensional but provided array has particle positions in {} dimensions",
                            System::DIM, arr.shape(1))
        );
    }
}

PYBIND11_MODULE(lv_mod, m) {
    ctiprd::binding::exportBaseTypes<System::dtype>(m);
    ctiprd::binding::exportSystem<System>(m, "LotkaVolterra");

    m.def("simulate", [] (std::size_t nSteps, float dt, int njobs,
            const np_array<System::dtype>& prey, const np_array<System::dtype>& predator,
            py::handle progressCallback) {
        check_shape(predator);
        check_shape(prey);

        System system {};

        ctiprd::binding::Trajectory<System> traj;

        auto pool = ctiprd::config::make_pool(njobs);
        auto integrator = ctiprd::cpu::integrator::EulerMaruyama{system, pool};
        for(auto i = 0; i < prey.shape(0); ++i) {
            integrator.particles()->addParticle({prey.at(i, 0), prey.at(i, 1)}, "prey");
        }
        for(auto i = 0; i < predator.shape(0); ++i) {
            integrator.particles()->addParticle({predator.at(i, 0), predator.at(i, 1)}, "predator");
        }

        {
            py::gil_scoped_release release;
            for (std::size_t step = 0; step < nSteps; ++step) {
                if (step % 20 == 0) {
                    traj.record(step, integrator, pool);
                }

                integrator.step(dt);

                if(step % 5 == 0) {
                    py::gil_scoped_acquire acquire;
                    if (PyErr_CheckSignals() != 0) {
                        throw py::error_already_set();
                    }

                    progressCallback(step);
                }
            }
        }

        pool->stop();
        return traj;
    });
}
