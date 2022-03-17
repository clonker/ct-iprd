#pragma once

#include <memory>

#include <ctiprd/util/distribution_utils.h>
#include <ctiprd/ParticleCollection.h>
#include <ctiprd/NeighborList.h>

namespace ctiprd::integrator {

namespace detail {
template <typename dtype>
auto &normalDistribution() {
    static thread_local std::normal_distribution<dtype> distribution {};
    return distribution;
}
}

template<int DIM, typename dtype, typename ExternalPotentials, typename PairPotentials,
         typename Pool = config::ThreadPool, typename Generator = std::mt19937>
class EulerMaruyama {
public:

    using Particles = ParticleCollection<DIM, dtype>;
    using NeighborList = nl::NeighborList<DIM, false, dtype, Pool>;
    static constexpr const char* name = "EulerMaruyama";

    static constexpr int nExternalPotentials = std::tuple_size_v<ExternalPotentials>;
    static constexpr int nPairPotentials = std::tuple_size_v<PairPotentials>;

    explicit EulerMaruyama(config::PoolPtr<Pool> pool) :
        particles_(std::make_shared<Particles>("A")), pool_(pool), externalPotentials_(), pairPotentials_(),
        neighborList_({5., 5.}, .5, pool)
    {
    }

    std::shared_ptr<Particles> particles() const {
        return particles_;
    }

    void step(double h) {
        neighborList_.update(particles_);
        auto worker = [
                h,
                &pot = externalPotentials_,
                &potPair = pairPotentials_,
                &nl = neighborList_,
                &data = *particles_
                        ]
                (auto id, typename Particles::Position &pos, typename Particles::Force &force) {

            force = {};
            std::apply([&pos, &force](auto &&... args) { ((force += args.force(pos)), ...); }, pot);

            nl.template forEachNeighbor(id, data, [&pos, &force, &potPair](auto neighborId, const auto &neighborPos, const auto &neighborForce) {
                std::apply([&pos, &force, &neighborPos](auto &&... args) { ((force += args.force(pos, neighborPos)), ...); }, potPair);
            });

            auto diffusionConstant = static_cast<dtype>(1);
            auto kbt = static_cast<dtype>(1);
            auto deterministicDisplacement = force * diffusionConstant * h / kbt;
            auto randomDisplacement = noise() * std::sqrt(2 * diffusionConstant * h);
            pos += deterministicDisplacement + randomDisplacement;
        };

        particles_->template forEachParticle(worker, pool_);
        pool_->waitForTasks();
    }

private:

    static typename Particles::Position noise() {
        typename Particles::Position out;
        std::generate(begin(out.data), end(out.data), [](){
            return detail::normalDistribution<dtype>()(rnd::staticThreadLocalGenerator<Generator>());
        });
        return out;
    }

    std::shared_ptr<Particles> particles_;
    config::PoolPtr<Pool> pool_;

    ExternalPotentials externalPotentials_;
    PairPotentials pairPotentials_;
    NeighborList neighborList_;
};
}
