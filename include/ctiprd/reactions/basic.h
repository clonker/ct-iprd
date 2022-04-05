/**
 *
 *
 * @file basic.h
 * @brief 
 * @author clonker
 * @date 3/27/22
 */
#pragma once

#include <ctiprd/util/distribution_utils.h>

namespace ctiprd::reactions {

template<typename dtype, typename Reactions>
dtype reactionRadius(const Reactions &reactions) {
    dtype result {};
    std::apply([&result](auto &&... args) { ((result = std::max(result, args.reactionRadius)), ...); }, reactions);
    return result;
}

namespace tags{

struct decay {};
struct conversion {};
struct fission {};
struct fusion {};
struct catalytic {};

}

namespace doi {

template<typename dtype>
struct Decay {
    using type = tags::decay;
    constexpr static std::size_t N_EDUCTS = 1;
    constexpr static std::size_t N_PRODUCTS = 0;

    std::size_t eductType;
    dtype rate;
};

template<typename dtype>
struct Conversion{
    using type = tags::conversion;
    constexpr static std::size_t N_EDUCTS = 1;
    constexpr static std::size_t N_PRODUCTS = 1;

    std::size_t eductType;
    std::size_t productType;
    dtype rate;
};

template<typename dtype>
struct Fission {
    using type = tags::fission;
    constexpr static std::size_t N_EDUCTS = 1;
    constexpr static std::size_t N_PRODUCTS = 2;

    std::size_t eductType;
    std::size_t productType1, productType2;
    dtype productDistance;
    dtype rate;
};

template<typename dtype>
struct Fusion {
    using type = tags::fusion;
    constexpr static std::size_t N_EDUCTS = 2;
    constexpr static std::size_t N_PRODUCTS = 1;

    std::size_t eductType1;
    std::size_t eductType2;
    std::size_t productType;
    dtype reactionRadius;
    dtype rate;
    dtype w1 {.5};
    dtype w2 {.5};
};

template<typename dtype>
struct Catalysis {
    using type = tags::catalytic;
    constexpr static std::size_t N_EDUCTS = 2;
    constexpr static std::size_t N_PRODUCTS = 2;

    std::size_t catalyst;
    std::size_t eductType;
    std::size_t productType;
    dtype reactionRadius;
    dtype rate;
};

}
}
