#pragma once

#include <memory>

#include <ctiprd/util/distribution_utils.h>
#include <ctiprd/ParticleCollection.h>
#include <ctiprd/NeighborList.h>
#include <ctiprd/potentials/util.h>
#include <ctiprd/potentials/ForceField.h>
#include <ctiprd/reactions/UncontrolledApproximation.h>

namespace ctiprd::integrator {

namespace detail {
template<typename dtype>
auto &normalDistribution() {
    static thread_local std::normal_distribution<dtype> distribution{};
    return distribution;
}

}

template<typename System, typename Pool = config::ThreadPool, typename Generator = std::mt19937,
         typename ForceField = potentials::ForceField<System>, typename Reactions = reactions::UncontrolledApproximation<System, Generator>>
class EulerMaruyama {
public:

    static constexpr std::size_t DIM = System::DIM;
    using dtype = typename System::dtype;

    using Particles = ParticleCollection<System>;
    static constexpr const char *name = "EulerMaruyama";

    explicit EulerMaruyama(const System &system, config::PoolPtr<Pool> pool) :
            particles_(std::make_shared<Particles>()), pool_(pool), system(system),
            forceField(std::make_unique<ForceField>()),
            reactions(std::make_unique<Reactions>()) {
    }

    std::shared_ptr<Particles> particles() const {
        return particles_;
    }

    void step(double h) {

        forceField->forces(particles_, pool_);

        auto worker = [h, &data = *particles_ ]
                (auto id, typename Particles::Position &pos, const typename Particles::ParticleType &type,
                 typename Particles::Force &force) {

            const auto &diffusionConstant = System::types[type].diffusionConstant;
            auto deterministicDisplacement = force * diffusionConstant * h;
            auto randomDisplacement = noise() * std::sqrt(2 * diffusionConstant * h);
            pos += deterministicDisplacement + randomDisplacement;
        };

        particles_->forEachParticle(worker, pool_);
        pool_->waitForTasks();

        reactions->reactions(system, h, particles_, pool_);
    }



private:

    static typename Particles::Position noise() {
        typename Particles::Position out;
        std::generate(begin(out.data), end(out.data), []() {
            return detail::normalDistribution<dtype>()(rnd::staticThreadLocalGenerator<Generator>());
        });
        return out;
    }

    std::shared_ptr<Particles> particles_;
    std::unique_ptr<ForceField> forceField;
    std::unique_ptr<Reactions> reactions;
    config::PoolPtr<Pool> pool_;
    System system;
};
}
