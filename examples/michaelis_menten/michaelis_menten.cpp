#include <cstddef>
#include <tuple>
#include <optional>

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#include <ctiprd/integrators/EulerMaruyama.h>
#include <ctiprd/potentials/external.h>
#include <ctiprd/ParticleTypes.h>
#include <ctiprd/reactions/basic.h>
#include <ctiprd/systems/util.h>
#include <ctiprd/binding/system_bindings.h>

namespace py = pybind11;

template<typename T>
struct MichaelisMenten {
    using dtype = T;
    static constexpr std::size_t DIM = 3;
    static constexpr std::array<T, DIM> boxSize{.3, .3, .3};
    static constexpr bool periodic = true;
    static constexpr ctiprd::ParticleTypes<dtype, 4> types{{
        { .name = "E", .diffusionConstant = 10.},
        { .name = "S", .diffusionConstant = 10.},
        { .name = "ES", .diffusionConstant = 10.},
        { .name = "P", .diffusionConstant = 10.}
        }};
    static constexpr std::size_t eId = ctiprd::systems::particleTypeId<types>("E");
    static constexpr std::size_t sId = ctiprd::systems::particleTypeId<types>("S");
    static constexpr std::size_t esId = ctiprd::systems::particleTypeId<types>("ES");
    static constexpr std::size_t pId = ctiprd::systems::particleTypeId<types>("P");

    MichaelisMenten() {
        {
            std::get<0>(reactionsO1) = ctiprd::reactions::doi::Fission<T>{
                    .eductType = esId,
                    .productType1 = eId,
                    .productType2 = sId,
                    .productDistance = .03,
                    .rate = 1.
            };
            std::get<0>(reactionsO1) = ctiprd::reactions::doi::Fission<T>{
                    .eductType = esId,
                    .productType1 = eId,
                    .productType2 = pId,
                    .productDistance = .03,
                    .rate = 1.
            };
        }
        {
            std::get<0>(reactionsO2) = ctiprd::reactions::doi::Fusion<T>{
                    .eductType1 = eId,
                    .eductType2 = sId,
                    .productType = esId,
                    .reactionRadius = 0.03,
                    .rate = 86.78638438
            };

        }
    }


    using ExternalPotentials = std::tuple<>;
    using PairPotentials = std::tuple<>;

    using ReactionsO1 = std::tuple<ctiprd::reactions::doi::Fission<T>, ctiprd::reactions::doi::Fission<T>>;
    using ReactionsO2 = std::tuple<ctiprd::reactions::doi::Fusion<T>>;

    ReactionsO1 reactionsO1{};
    ReactionsO2 reactionsO2{};

    ExternalPotentials externalPotentials{};
    PairPotentials pairPotentials{};
};

using System = MichaelisMenten<float>;

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
                integrator.step(dt);

                nE.emplace_back();
                nS.emplace_back();
                nES.emplace_back();
                nP.emplace_back();

                integrator.particles()->forEachParticle([&m, &nE, &nS, &nES, &nP](auto id, const auto &pos, const auto &type,
                                                           const auto &force) {
                    std::size_t nEl {};
                    std::size_t nSl {};
                    std::size_t nESl {};
                    std::size_t nPl {};

                    if(type == System::eId) {
                        ++nEl;
                    } else if(type == System::sId) {
                        ++nSl;
                    } else if(type == System::esId) {
                        ++nESl;
                    } else if(type == System::pId) {
                        ++nPl;
                    }

                    {
                        std::scoped_lock lock {m};
                        nE.back() += nEl;
                        nS.back() += nSl;
                        nES.back() += nESl;
                        nP.back() += nPl;
                    }
                }, pool);
                pool->waitForTasks();

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
