#pragma once

#include <ctiprd/util/distribution_utils.h>
#include <memory>
#include <ctiprd/ParticleCollection.h>

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
    static constexpr const char* name = "EulerMaruyama";

    static constexpr int nExternalPotentials = std::tuple_size_v<ExternalPotentials>;
    static constexpr int nPairPotentials = std::tuple_size_v<PairPotentials>;

    explicit EulerMaruyama(config::PoolPtr<Pool> pool) : particles_(std::make_shared<Particles>("A")), pool_(pool) { }

    std::shared_ptr<Particles> particles() const {
        return particles_;
    }

    void forces() {
        auto workerExternal = [](auto index, const typename Particles::MaybePosition &pos, typename Particles::Force &force) {
            force = {};
            if (pos) {
                std::apply([&pos, &force](auto &&... args) { ((force += args.force(*pos)), ...); }, ExternalPotentials{});
            }
        };
        auto futures = particles_->template forEachParticle(workerExternal, pool_);
        std::for_each(begin(futures), end(futures), [](auto &f) { f.wait(); });
    }

    void step(double h) {
        auto worker = [h](auto index, typename Particles::MaybePosition &pos, const typename Particles::Force &force) {
            if (pos) {
                // if the particle exists
                auto diffusionConstant = static_cast<dtype>(1);
                auto kbt = static_cast<dtype>(1);
                auto deterministicDisplacement = force * diffusionConstant * h / kbt;
                auto randomDisplacement = noise() * std::sqrt(2 * diffusionConstant * h);
                *pos += deterministicDisplacement + randomDisplacement;
            }
        };

        auto futures = particles_->template forEachParticle(worker, pool_);
        std::for_each(begin(futures), end(futures), [](auto &f) { f.wait(); });
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
};
}
