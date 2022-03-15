//
// Created by mho on 3/15/22.
//

#include <cmath>

#include "ctiprd/ParticleCollection.h"

static constexpr int DIM = 2;
using dtype = float;

using Collection = ctiprd::ParticleCollection<DIM, dtype>;

struct DoubleWell2D {
    [[nodiscard]] static constexpr dtype energy(const Collection::Position &x) {
        return (x[0] * x[0] - static_cast<dtype>(1)) * (x[0] * x[0] - static_cast<dtype>(1)) + x[1] * x[1];
    }

    [[nodiscard]] static constexpr Collection::Force f(const Collection::Position &x) {
        return {-4 * x[0] * x[0] * x[0] + 4 * x[0], -2 * x[1]};
    }

    dtype mass{1.};
    dtype damping{1.};
    dtype kT{1.};

    dtype sigma {std::sqrt(static_cast<dtype>(0.5) * kT / (mass * damping))};
};

int main() {
    auto pool = ctiprd::config::make_pool(5);
    return 0;
}
