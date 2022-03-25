//
// Created by mho on 3/17/22.
//

#pragma once

#include <cstddef>
#include <cmath>

#include <ctiprd/vec.h>

namespace ctiprd::potential::pair {

template<typename dtype>
struct HarmonicRepulsion {

    static constexpr std::size_t DIM = 2;
    using State = Vec<dtype, DIM>;

    dtype energy(const State &x1, const State &x2) const {
        auto dSquared = (x2 - x1).normSquared();
        if (dSquared < cutoff * cutoff) {
            auto d = std::sqrt(dSquared);
            d -= cutoff;
            return static_cast<dtype>(0.5) * d * d * forceConstant;
        }
        return 0;
    }

    State force(const State &x1, const State &x2) const {
        auto xij = x2 - x1;
        auto dSquared = xij.normSquared();
        if (dSquared < cutoff * cutoff && dSquared > 0) {
            auto d = std::sqrt(dSquared);
            return (forceConstant * (d - cutoff)) / d * xij;
        } else {
            return {};
        }
    }

    dtype cutoff{1.};
    dtype forceConstant{1.};
};
}
