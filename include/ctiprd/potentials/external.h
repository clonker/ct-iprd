//
// Created by mho on 3/16/22.
//
#pragma once

#include <cstddef>
#include <cmath>

#include <ctiprd/geometry/box.h>

namespace ctiprd::potentials::external {

template<typename Geometry, std::size_t dim, bool inclusion, bool allTypes, std::size_t typeId = 0>
struct HarmonicGeometry {
    using dtype = typename Geometry::dtype;
    constexpr static std::size_t particleType = typeId;
    static constexpr std::size_t DIM = dim;

    template<typename State, typename ParticleType>
    [[nodiscard]] constexpr dtype energy(const State &x, const ParticleType &type) const {
        if (allTypes || type == particleType) {
            auto shortestDiff = geometry.template smallestDifference<inclusion>(x);
            return 0.5 * k * shortestDiff.normSquared();
        }
        return {};
    }

    template<typename State, typename ParticleType>
    [[nodiscard]] constexpr State force(const State &x, const ParticleType &type) const {
        if (allTypes || type == particleType) {
            return -1 * k * geometry.template smallestDifference<inclusion>(x);
        }
        return {};
    }

    dtype k {1.};
    Geometry geometry {};
};

template<typename dtype, std::size_t dim, bool allTypes, std::size_t typeId = 0>
using BoxInclusion = HarmonicGeometry<geometry::Box<dim, dtype>, dim, true, allTypes, typeId>;
template<typename dtype, std::size_t dim, bool allTypes, std::size_t typeId = 0>
using BoxExclusion = HarmonicGeometry<geometry::Box<dim, dtype>, dim, false, allTypes, typeId>;

template<typename dtype, bool allTypes, std::size_t typeId = 0>
struct DoubleWell {

    constexpr static std::size_t particleType = typeId;

    static constexpr std::size_t DIM = 2;

    template<typename State, typename ParticleType>
    [[nodiscard]] constexpr dtype energy(const State &x, const ParticleType &type) const {
        if (allTypes || type == particleType) {
            return k*(x[0] * x[0] - 1.) * (x[0] * x[0] - 1.) + k*x[1] * x[1];
        }
        return {};
    }

    template<typename State, typename ParticleType>
    [[nodiscard]] constexpr State force(const State &x, const ParticleType &type) const {
        if (allTypes || type == particleType) {
            return {{-4 * k * x[0] * x[0] * x[0] + 4 * k * x[0], -2 * k * x[1]}};
        }
        return {};
    }

    dtype k {1.};
};

}
