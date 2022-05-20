//
// Created by mho on 5/20/22.
//

#include <cstddef>
#include <tuple>
#include <optional>

#include <ctiprd/vec.h>
#include <ctiprd/potentials/external.h>
#include <ctiprd/potentials/interaction.h>
#include <ctiprd/ParticleTypes.h>
#include <ctiprd/reactions/doi.h>
#include <ctiprd/systems/util.h>
#include <ctiprd/util/rates.h>

#include <ctiprd/config.h>
#include <ctiprd/binding/system_bindings.h>
#include <ctiprd/progressbar.hpp>

#include <ctiprd/cpu/integrators/EulerMaruyama.h>

template<typename T>
struct Conf {
    static constexpr T diffPrey = 0.01;
    static constexpr T diffPred = 0.01;

    static constexpr T alpha = 2.;  // birth: prey -> prey + prey
    static constexpr T alphaDistance = .2;  // birth distance
    static constexpr T beta = 0.05;  // eat: prey + pred -> pred + pred
    static constexpr T betaRadius = 0.25;
    static constexpr T betaMic = 7.670679846561291;
    static constexpr T gamma = 1.5;

    static constexpr T friction = 0.01;
    static constexpr T frictionRadius = 0.2;
    static constexpr T frictionMic = 0.39155565024786165;
};

template<typename T>
struct LotkaVolterra2d {
    using Cfg = Conf<T>;
    using dtype = T;
    static constexpr std::size_t DIM = 2;
    static constexpr std::array<T, DIM> boxSize{10., 50.};
    static constexpr bool periodic = true;
    static constexpr T kBT = 2.43614;
    static constexpr ctiprd::ParticleTypes<dtype, 3> types{{
                                                                   {
                                                                           .name = "predator",
                                                                           .diffusionConstant = Cfg::diffPred
                                                                   },
                                                                   {
                                                                           .name = "prey",
                                                                           .diffusionConstant = Cfg::diffPrey
                                                                   },
                                                                   {
                                                                       .name = "barrier",
                                                                       .diffusionConstant = 0
                                                                   }
                                                           }};
    static constexpr std::size_t preyId = ctiprd::systems::particleTypeId<types>("prey");
    static constexpr std::size_t predatorId = ctiprd::systems::particleTypeId<types>("predator");
    static constexpr std::size_t wallId = ctiprd::systems::particleTypeId<types>("barrier");

    LotkaVolterra2d() {
        {
            auto &[birth, death] = reactionsO1;
            birth = ctiprd::reactions::doi::Fission<T>{
                    .eductType = preyId,
                    .productType1 = preyId,
                    .productType2 = preyId,
                    .productDistance = Cfg::alphaDistance,
                    .rate = Cfg::alpha
            };
            death = ctiprd::reactions::doi::Decay<T>{
                    .eductType = predatorId,
                    .rate = Cfg::gamma
            };
        }
        {
            auto &[preySocialFriction, predatorSocialFriction, predEatsPrey] = reactionsO2;
            preySocialFriction = ctiprd::reactions::doi::Fusion<T>{
                    .eductType1 = preyId,
                    .eductType2 = preyId,
                    .productType = preyId,
                    .reactionRadius = Cfg::frictionRadius,
                    .rate = Cfg::frictionMic
            };
            predatorSocialFriction = ctiprd::reactions::doi::Fusion<T>{
                    .eductType1 = predatorId,
                    .eductType2 = predatorId,
                    .productType = predatorId,
                    .reactionRadius = Cfg::frictionRadius,
                    .rate = Cfg::frictionMic
            };
            predEatsPrey = ctiprd::reactions::doi::Catalysis<T>{
                    .catalyst = predatorId,
                    .eductType = preyId,
                    .productType = predatorId,
                    .reactionRadius = Cfg::betaRadius,
                    .rate = Cfg::betaMic
            };
        }


        {
            static constexpr dtype eps = 1e-5;

            auto kmac = ctiprd::rates::macroscopicRate(Cfg::betaMic, types[preyId].diffusionConstant,
                                                       types[predatorId].diffusionConstant, Cfg::betaRadius);
            if (std::abs(kmac - Cfg::beta) > eps) {
                throw std::runtime_error(fmt::format("kmac = {}, beta = {}", kmac, Cfg::beta));
            }
        }
        {
            static constexpr dtype eps = 1e-5;

            auto kmac = ctiprd::rates::macroscopicRate(Cfg::frictionMic, types[preyId].diffusionConstant,
                                                       types[predatorId].diffusionConstant, Cfg::frictionRadius);
            if (std::abs(kmac - Cfg::friction) > eps) {
                throw std::runtime_error(fmt::format("kmac = {}, friction = {}", kmac, Cfg::friction));
            }
        }
        {
            auto &box = std::get<0>(externalPotentials);
            box.geometry.v0 = {-4., -24.};
            box.geometry.v1 = {4., 24.};
            box.k = 150.;
        }
        {
            auto &[preyWall, predWall] = pairPotentials;

            preyWall.cutoff = .3;
            preyWall.forceConstant = 50.;
            preyWall.particleType1 = wallId;
            preyWall.particleType2 = preyId;

            predWall.cutoff = .3;
            predWall.forceConstant = 50.;
            predWall.particleType1 = wallId;
            predWall.particleType2 = predatorId;
        }
    }


    using ExternalPotentials = std::tuple<ctiprd::potentials::external::BoxInclusion<dtype, DIM>>;
    using PairPotentials = std::tuple<ctiprd::potentials::pair::HarmonicRepulsion<dtype>,
                                      ctiprd::potentials::pair::HarmonicRepulsion<dtype>>;

    using ReactionsO1 = std::tuple<
            ctiprd::reactions::doi::Fission<T>, // prey birth
            ctiprd::reactions::doi::Decay<T> // predator death
    >;
    using ReactionsO2 = std::tuple<
            ctiprd::reactions::doi::Fusion<T>, // prey social friction
            ctiprd::reactions::doi::Fusion<T>, // predator social friction
            ctiprd::reactions::doi::Catalysis<T> // predator eats prey
    >;

    ReactionsO1 reactionsO1{};
    ReactionsO2 reactionsO2{};

    ExternalPotentials externalPotentials{};
    PairPotentials pairPotentials{};
};

using System = LotkaVolterra2d<float>;
namespace py = pybind11;

template<typename T>
using np_array = ctiprd::binding::np_array<T>;

void check_shape(const np_array <System::dtype> &arr) {
    if (arr.ndim() != 2) {
        throw std::runtime_error(
                fmt::format("Provided particle array needs to be two-dimensional but was {}-dimensional.", arr.ndim())
        );
    }
    if (arr.shape(1) != System::DIM) {
        throw std::runtime_error(
                fmt::format("System is {}-dimensional but provided array has particle positions in {} dimensions",
                            System::DIM, arr.shape(1))
        );
    }
}

PYBIND11_MODULE(lv2d_mod, m) {
    ctiprd::binding::exportBaseTypes<System::dtype>(m);
    ctiprd::binding::exportSystem<System>(m, "LotkaVolterra");

    m.def("simulate", [](std::size_t nSteps, float dt, int njobs,
                         const np_array <System::dtype> &prey, const np_array <System::dtype> &predator,
                         const np_array<System::dtype> &walls,
                         py::handle progressCallback) {
        check_shape(predator);
        check_shape(prey);
        check_shape(walls);

        System system{};

        ctiprd::binding::Trajectory<System> traj;

        auto pool = ctiprd::config::make_pool(njobs);
        auto integrator = ctiprd::cpu::integrator::EulerMaruyama{system, pool};
        for (auto i = 0; i < prey.shape(0); ++i) {
            integrator.particles()->addParticle({prey.at(i, 0), prey.at(i, 1)}, "prey");
        }
        for (auto i = 0; i < predator.shape(0); ++i) {
            integrator.particles()->addParticle({predator.at(i, 0), predator.at(i, 1)}, "predator");
        }

        for (auto i = 0; i < walls.shape(0); ++i) {
            integrator.particles()->addParticle({walls.at(i, 0), walls.at(i, 1)}, "barrier");
        }

        {
            py::gil_scoped_release release;
            for (std::size_t step = 0; step < nSteps; ++step) {
                if (step % 20 == 0) {
                    traj.record(step, integrator, pool);
                }

                integrator.step(dt);

                if (step % 5 == 0) {
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
