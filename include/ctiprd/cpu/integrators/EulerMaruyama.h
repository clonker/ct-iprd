#pragma oncecpu

#include <memory>

#include <ctiprd/potentials/util.h>
#include <ctiprd/util/distribution_utils.h>

#include <ctiprd/cpu/ParticleCollection.h>
#include <ctiprd/cpu/NeighborList.h>
#include <ctiprd/cpu/ForceField.h>
#include <ctiprd/cpu/UncontrolledApproximation.h>
#include <ctiprd/util/pbc.h>

namespace ctiprd::cpu::integrator {

namespace detail {
template<typename dtype>
auto &normalDistribution() {
    static thread_local std::normal_distribution<dtype> distribution{};
    return distribution;
}

}

template<typename System, typename Pool = config::ThreadPool, typename Generator = std::mt19937,
         typename ForceField = potentials::ForceField<System>,
         typename ParticleCollection = ParticleCollection<System, particles::positions, particles::forces>,
         typename Reactions = UncontrolledApproximation<ParticleCollection, System, Generator>>
class EulerMaruyama {
public:
    static constexpr std::size_t DIM = System::DIM;
    using Particles = ParticleCollection;
    using dtype = typename System::dtype;

    static constexpr const char *name = "EulerMaruyama";

    explicit EulerMaruyama(const System &system, config::PoolPtr<Pool> pool) :
            particles_(std::make_shared<Particles>()), pool_(pool), system(system),
            forceField(std::make_unique<ForceField>()),
            reactions(std::make_unique<Reactions>(system)) {
        forceField->setExternalPotentials(system.externalPotentials);
        forceField->setPairPotentials(system.pairPotentials);
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
            auto deterministicDisplacement = force * diffusionConstant * h / System::kBT;
            auto randomDisplacement = noise() * std::sqrt(2 * diffusionConstant * h);
            pos += deterministicDisplacement + randomDisplacement;

            util::pbc::wrapPBC<System>(pos);
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
