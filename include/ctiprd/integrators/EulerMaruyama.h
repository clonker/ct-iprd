#pragma once

#include <memory>

#include <ctiprd/util/distribution_utils.h>
#include <ctiprd/ParticleCollection.h>
#include <ctiprd/NeighborList.h>
#include <ctiprd/potentials/util.h>

namespace ctiprd::integrator {

namespace detail {
template<typename dtype>
auto &normalDistribution() {
    static thread_local std::normal_distribution<dtype> distribution{};
    return distribution;
}

}

template<typename System,
        typename Pool = config::ThreadPool, typename Generator = std::mt19937>
class EulerMaruyama {
public:

    static constexpr std::size_t DIM = System::DIM;
    using dtype = typename System::dtype;
    using ExternalPotentials = typename System::ExternalPotentials;
    using PairPotentials = typename System::PairPotentials;

    using Particles = ParticleCollection<System>;
    using NeighborList = nl::NeighborList<DIM, false, dtype, Pool>;
    static constexpr const char *name = "EulerMaruyama";

    static constexpr int nExternalPotentials = std::tuple_size_v<ExternalPotentials>;
    static constexpr int nPairPotentials = std::tuple_size_v<PairPotentials>;

    explicit EulerMaruyama(config::PoolPtr<Pool> pool) :
            particles_(std::make_shared<Particles>()), pool_(pool), externalPotentials_(), pairPotentials_() {
    }

    std::shared_ptr<Particles> particles() const {
        return particles_;
    }

    void step(double h) {
        if (!neighborList_) {
            auto cutoff = potentials::cutoff<dtype>(pairPotentials_);
            neighborList_ = std::make_unique<NeighborList>(System::boxSize, cutoff, pool_);
        }
        neighborList_->update(particles_);
        forces();

        auto worker = [
                h,
                &pot = externalPotentials_,
                &potPair = pairPotentials_,
                nl = neighborList_.get(),
                &data = *particles_
        ]
                (auto id, typename Particles::Position &pos, const typename Particles::ParticleType &type,
                 typename Particles::Force &force) {

            const auto &diffusionConstant = System::types[type].diffusionConstant;
            auto deterministicDisplacement = force * diffusionConstant * h;
            auto randomDisplacement = noise() * std::sqrt(2 * diffusionConstant * h);
            pos += deterministicDisplacement + randomDisplacement;
        };

        particles_->template forEachParticle(worker, pool_);
        pool_->waitForTasks();
    }

    void setExternalPotentials(const ExternalPotentials &potentials) {
        externalPotentials_ = potentials;
    }

    void setPairPotentials(const PairPotentials &potentials) {
        pairPotentials_ = potentials;
        neighborList_.reset();
    }

private:

    void forces() {
        auto worker = [
                &pot = externalPotentials_,
                &potPair = pairPotentials_,
                nl = neighborList_.get(),
                &data = *particles_
        ]
                (auto id, typename Particles::Position &pos, const typename Particles::ParticleType &type,
                 typename Particles::Force &force) {

            std::fill(begin(force.data), end(force.data), static_cast<dtype>(0));
            std::apply([&pos, &type, &force](auto &&... args) {
                ((force += args.force(pos, type)), ...);
            }, pot);

            nl->template forEachNeighbor(id, data, [&pos, &type, &force, &potPair](auto neighborId, const auto &neighborPos, const auto &neighborType,
                                                                            const auto &neighborForce) {
                std::apply([&pos, &type, &force, &neighborPos, &neighborType](auto &&... args) {
                    ((force += args.force(pos, type, neighborPos, neighborType)), ...);
                }, potPair);
            });
        };
        particles_->template forEachParticle(worker, pool_);
        pool_->waitForTasks();
    }

    static typename Particles::Position noise() {
        typename Particles::Position out;
        std::generate(begin(out.data), end(out.data), []() {
            return detail::normalDistribution<dtype>()(rnd::staticThreadLocalGenerator<Generator>());
        });
        return out;
    }

    std::shared_ptr<Particles> particles_;
    config::PoolPtr<Pool> pool_;

    ExternalPotentials externalPotentials_;
    PairPotentials pairPotentials_;
    std::unique_ptr<NeighborList> neighborList_;
};
}
