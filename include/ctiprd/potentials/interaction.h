//
// Created by mho on 3/17/22.
//

#pragma once

#include <cstddef>
#include <cmath>

#include <ctiprd/util/pbc.h>
#include <ctiprd/potentials/util.h>

namespace ctiprd::potentials::pair {

template<typename dtype>
struct HarmonicRepulsion {

    template<typename ParticleType>
    [[nodiscard]] bool supportsTypes(const ParticleType &type1, const ParticleType &type2) const {
        return particleType1 == -1 || (type1 == particleType1 && type2 == particleType2) ||
                (type1 == particleType2 && type2 == particleType1);
    }

    dtype cutoff{1.};
    dtype forceConstant{1.};

    int particleType1 {-1};
    int particleType2 {-1};
};
}
