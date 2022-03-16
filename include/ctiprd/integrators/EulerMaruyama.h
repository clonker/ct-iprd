#pragma once

#include <ctiprd/util/distribution_utils.h>
#include <memory>
#include <ctiprd/ParticleCollection.h>

namespace ctiprd::integrator {

template<int DIM, typename dtype, typename ExternalPotentials, typename PairPotentials, typename Generator = std::mt19937>
class EulerMaruyama {
public:

    using Particles = ParticleCollection<DIM, dtype>;

    static constexpr const char* name = "EulerMaruyama";

    explicit EulerMaruyama(std::shared_ptr<> particles) : particles(particles) {
    }

    /*template<typename F, typename Sigma>
    void eval(F &&f, const Sigma& sigma, double h, std::size_t nSteps, double t0, const State &y0) {
        auto y = y0;
        auto t = t0;
        auto sqrth = std::sqrt(h);
        for (std::size_t i = 0; i < nSteps; ++i) {
            y = step(f, h, sqrth, t, y, sigma);
            t += h;
        }

        return y;
    }

    template<typename F, typename Sigma>
    State step(F &&f, double h, double sqrth, double t, const State &y, const Sigma &sigma) {
        auto mu = f(t, y);
        auto w = noise(); // evaluate Wiener processes
        auto out = y;

        for (std::size_t j = 0; j < DIM; ++j) {
            out[j] += h * mu[j];

            for (std::size_t k = 0; k < DIM; ++k) {
                out[j] += sigma[j][k] * sqrth * w[k];
            }
        }
        return out;
    }*/

private:

    /*State noise() {
        State out;
        std::generate(begin(out), end(out), [this](){
            return distribution(generator);
        });
        return out;
    }*/

    Generator generator;
    std::normal_distribution<Value> distribution;

    std::shared_ptr<ParticleCollection> particles;
};
}
