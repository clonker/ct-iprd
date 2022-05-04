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
    using Info = systems::SystemInfo<System>;
    using Particles = ParticleCollection;

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

    void step(double stepSize) {

        if(prevIntegrationStep != stepSize) {
            prevIntegrationStep = stepSize;

            for(auto i = 0; i < Info::nTypes; ++i) {
                randomDisplacementPrefactors[i] = std::sqrt(2 * Info::diffusionConstantOf(i) * stepSize);
                deterministicDisplacementPrefactors[i] = Info::diffusionConstantOf(i) * stepSize / System::kBT;
            }
        }

        if constexpr(Info::hasForces()) {
            forceField->forces(particles_, pool_);
        }

        const auto worker = [this]
                (const auto &, typename Particles::Position &pos, const typename Particles::ParticleType &type,
                 const typename Particles::Force &force) {
            pos += force * deterministicDisplacementPrefactors[type] + noise() * randomDisplacementPrefactors[type];
            util::pbc::wrapPBC<System>(pos);
        };

        auto futures = particles_->forEachParticle(worker, pool_);
        for(auto &future : futures) {
            future.wait();
        }

        if constexpr(Info::hasReactions()) {
            reactions->reactions(stepSize, particles_, pool_);

            if constexpr(Info::periodic) {
                const auto fut = particles_->forEachParticle([](const auto &, auto &pos, const auto &, const auto &) {
                    util::pbc::wrapPBC<System>(pos);
                }, pool_);
                for(auto &future : fut) {
                    future.wait();
                }
            }
        }
    }

private:

    static typename Particles::Position noise() {
        typename Particles::Position out;
        std::generate(begin(out.data), end(out.data), []() {
            return detail::normalDistribution<typename Info::dtype>()(rnd::staticThreadLocalGenerator<Generator>());
        });
        return out;
    }

    std::shared_ptr<Particles> particles_;
    std::unique_ptr<ForceField> forceField;
    std::unique_ptr<Reactions> reactions;
    std::array<typename Info::dtype, Info::nTypes> randomDisplacementPrefactors {};
    std::array<typename Info::dtype, Info::nTypes> deterministicDisplacementPrefactors {};
    typename Info::dtype prevIntegrationStep {0};
    config::PoolPtr<Pool> pool_;
    System system;
};
}
