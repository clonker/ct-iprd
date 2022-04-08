//
// Created by mho on 3/17/22.
//

#pragma once

#include <cstddef>
#include <cmath>

#include <ctiprd/potentials/util.h>

namespace ctiprd::potentials::pair {

template<typename System, bool allTypes, std::size_t type1=0, std::size_t type2=0>
struct HarmonicRepulsion {
    using dtype = typename System::dtype;

    static constexpr std::size_t particleType1 = type1;
    static constexpr std::size_t particleType2 = type2;

    template<typename State, typename ParticleType>
    dtype energy(const State &x1, const ParticleType &t1, const State &x2, const ParticleType &t2) const {
        if(allTypes || potentialApplicable<HarmonicRepulsion>(t1, t2)) {
            auto dSquared = util::pbc::dSquared<System>(x1, x2);
            if (dSquared < cutoff * cutoff) {
                auto d = std::sqrt(dSquared);
                d -= cutoff;
                return static_cast<dtype>(0.5) * d * d * forceConstant;
            }
        }
        return 0;
    }

    template<typename State, typename ParticleType>
    State force(const State &x1, const ParticleType &t1, const State &x2, const ParticleType &t2) const {
        if (allTypes || potentialApplicable<HarmonicRepulsion>(t1, t2)) {
            auto xij = util::pbc::shortestDifference<System>(x1, x2);
            auto dSquared = xij.normSquared();
            if (dSquared < cutoff * cutoff && dSquared > 0) {
                auto d = std::sqrt(dSquared);
                return (forceConstant * (d - cutoff)) / d * xij;
            }
        }
        return {};
    }

    dtype cutoff{1.};
    dtype forceConstant{1.};
};
}
