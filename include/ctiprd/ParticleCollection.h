//
// Created by mho on 5/5/21.
//

#pragma once

#include <cstdint>
#include <vector>

namespace ctiprd {

using ParticleType = std::uint32_t;

template<int DIM, typename dtype>
class ParticleCollection {
public:
    static constexpr int dim = DIM;

    ParticleCollection() = default;

    ParticleCollection(std::vector<dtype> particles, std::size_t nParticles)
            : _nParticles(nParticles), _particles(std::move(particles)) {
        _forces.resize(_particles.size());
        _types.resize(_particles.size());
    }

    [[nodiscard]] std::size_t nParticles() const { return _nParticles; }

    dtype *positions() {
        return _particles.data();
    }

    [[nodiscard]] const dtype *positions() const {
        return _particles.data();
    }

    dtype *position(std::size_t id) {
        return _particles.data() + id * DIM;
    }

    [[nodiscard]] const dtype *position(std::size_t id) const {
        return _particles.data() + id * DIM;
    }

    std::vector<dtype> &forces() { return _forces; }

    [[nodiscard]] const std::vector<dtype> &forces() const { return _forces; }

    dtype *force(std::size_t id) {
        return _forces.data() + id * DIM;
    }

    [[nodiscard]] const dtype *force(std::size_t id) const {
        return _forces.data() + id * DIM;
    }

    auto *types() { return _types.data(); }

    const auto *types() const { return _types.data(); }

    auto type(std::size_t id) { return _types[id]; }

    auto type(std::size_t id) const { return _types[id]; }

    template<typename F, typename Pool>
    void forEachParticle(F op, Pool &pool) {
        auto n = nParticles();
        auto *forces = this->forces().data();
        auto *pos = positions();

        auto granularity = config::threadGranularity(pool);
        auto grainSize = n / granularity;

        auto loop = [&op, pos, forces](auto beginIx, auto endIx) {
            for(std::size_t i = beginIx; i < endIx; ++i) {
                op(i, pos + i * dim, forces + i * dim);
            }
        };
        std::vector<std::future<void>> futures;
        futures.reserve(granularity);
        std::size_t beginIx = 0;
        for (std::size_t i = 0; i < granularity-1; ++i) {
            auto endIx = beginIx + grainSize;
            if(beginIx != endIx) {
                futures.push_back(pool.push(loop, beginIx, endIx));
            }
            beginIx = endIx;
        }
        if(beginIx != n) {
            futures.push_back(pool.push(loop, beginIx, n));
        }
        std::for_each(begin(futures), end(futures), [](auto &future) {future.wait();});
    }

private:
    std::size_t _nParticles{};
    std::vector<dtype> _particles{};
    std::vector<dtype> _forces{};
    std::vector<ParticleType> _types{};
};
}
