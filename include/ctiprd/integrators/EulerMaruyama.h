#pragma once

#include <ctiprd/util/distribution_utils.h>

namespace ctiprd::integrator {

template<typename State, std::size_t DIM, typename Value = double, typename Generator = std::mt19937>
class EulerMaruyama {
public:

    static constexpr const char* name = "EulerMaruyama";

    explicit EulerMaruyama(std::int64_t seed = -1) {
        generator = seed < 0 ? rnd::randomlySeededGenerator<Generator>() : rnd::seededGenerator<Generator>(seed);
    }

    template<typename F, typename Sigma>
    State eval(F &&f, const Sigma& sigma, double h, std::size_t nSteps, double t0, const State &y0) {
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
    }

private:

    State noise() {
        State out;
        std::generate(begin(out), end(out), [this](){
            return distribution(generator);
        });
        return out;
    }

    Generator generator;
    std::normal_distribution<Value> distribution;
};
}
