//
// Created by mho on 3/16/22.
//
#pragma once

#include <cstddef>
#include <cmath>

#include <ctiprd/geometry/box.h>

namespace ctiprd::potentials::external {

template<typename Geometry, bool inc>
struct HarmonicGeometry {

    static constexpr bool inclusion = inc;

    template<typename ParticleType>
    [[nodiscard]] bool supportsType(const ParticleType &type) const {
        return particleType == -1 || type == particleType;
    }

    typename Geometry::dtype k {1.};
    Geometry geometry {};
    int particleType {-1};
};

template<typename dtype, std::size_t dim>
using BoxInclusion = HarmonicGeometry<geometry::Box<dim, dtype>, true>;
template<typename dtype, std::size_t dim>
using BoxExclusion = HarmonicGeometry<geometry::Box<dim, dtype>, false>;

template<typename dtype>
struct DoubleWell {

    static constexpr std::size_t DIM = 2;

    template<typename ParticleType>
    [[nodiscard]] bool supportsType(const ParticleType &type) const {
        return particleType == -1 || type == particleType;
    }

    dtype k {1.};
    int particleType {-1};
};

}
