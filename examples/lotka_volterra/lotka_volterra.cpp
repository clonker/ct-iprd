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
using System = ctiprd::systems::LotkaVolterra<float>;

void check_shape(const np_array<float> &arr) {
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
            const np_array<float>& prey, const np_array<float>& predator,
            py::handle progressCallback) {
        check_shape(predator);
        check_shape(prey);

        System system {};
        auto pool = ctiprd::config::make_pool(njobs);
        auto integrator = ctiprd::cpu::integrator::EulerMaruyama{system, pool};
        for(auto i = 0; i < prey.shape(0); ++i) {
            integrator.particles()->addParticle({prey.at(i, 0), prey.at(i, 1)}, "prey");
        }
        for(auto i = 0; i < predator.shape(0); ++i) {
            integrator.particles()->addParticle({predator.at(i, 0), predator.at(i, 1)}, "predator");
        }

        std::vector<std::tuple<np_array<float>, np_array<std::size_t>>> trajectory;
        {
            py::gil_scoped_release release;
            for (std::size_t t = 0; t < nSteps; ++t) {
                integrator.step(dt);
                std::vector<std::future<void>> futures {};

                if (t % 200 == 0) {
                    py::gil_scoped_acquire acquire;

                    auto nParticles = integrator.particles()->size();
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
                    auto &[traj, types] = trajectory.back();

                    auto* trajBegin = traj.mutable_data(0);
                    auto* trajEnd = trajBegin + traj.size();
                    std::fill(trajBegin, trajEnd, std::numeric_limits<System::dtype>::infinity());

                    auto* typesBegin = types.mutable_data(0);
                    auto* typesEnd = typesBegin + types.size();
                    std::fill(typesBegin, typesEnd, std::numeric_limits<std::size_t>::infinity());

                    futures = integrator.particles()->forEachParticle([&trajBegin, &typesBegin](auto id, const auto &pos, const auto& type, const auto &force) {
                        typesBegin[id] = type;
                        trajBegin[2*id] = pos[0];
                        trajBegin[2*id+1] = pos[1];
                    }, pool);

                    progressCallback(t);
                }

                if(t % 50 == 0) {
                    py::gil_scoped_acquire acquire;
                    if (PyErr_CheckSignals() != 0) {
                        throw py::error_already_set();
                    }
                }

                for(auto &future : futures) future.wait();
            }
        }

        return trajectory;
    });
}
