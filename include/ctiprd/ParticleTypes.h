#pragma once

#include <array>
#include <string_view>

namespace ctiprd {

template<typename dtype>
struct ParticleType {
    std::string_view name;
    dtype diffusionConstant;
};

template<typename dtype, std::size_t N>
using ParticleTypes = std::array<ParticleType<dtype>, N>;

}
