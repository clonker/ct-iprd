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
struct ReactionE1 {
    constexpr static std::size_t N_EDUCTS = 1;

    template<typename State, typename ParticleType>
    [[nodiscard]] bool shouldPerform(dtype tau, const State &state, const ParticleType &t) const {
        return t == eductType && rnd::uniform_real<dtype>() < 1 - std::exp(-rate * tau);
    }

    std::size_t eductType;
    dtype rate;
};

template<typename dtype>
struct ReactionE2 {
    constexpr static std::size_t N_EDUCTS = 2;

    template<typename State, typename ParticleType>
    [[nodiscard]] bool shouldPerform(dtype tau, const State &state, const ParticleType &t,
                                     const State &state2, const ParticleType &t2) const {
        if ((t == particleType1 && t2 == particleType2) || (t == particleType2 && t2 == particleType1)) {
            return (state - state2).normSquared() < reactionRadius * reactionRadius &&
                   rnd::uniform_real<dtype>() < 1 - std::exp(-rate * tau);
        }
        return false;
    }

    std::size_t particleType1;
    std::size_t particleType2;
    dtype rate;
    dtype reactionRadius;
};

template<typename dtype>
struct Decay : public ReactionE1<dtype> {
    using type = tags::decay;
    constexpr static std::size_t N_PRODUCTS = 0;
};

template<typename dtype>
struct Conversion : public ReactionE1<dtype> {
    using type = tags::conversion;
    constexpr static std::size_t N_PRODUCTS = 1;

    std::size_t productType;
};

template<typename dtype>
struct Fission : public ReactionE1<dtype> {
    using type = tags::fission;
    constexpr static std::size_t N_PRODUCTS = 2;

    dtype productDistance;
};

template<typename dtype>
struct Fusion : public ReactionE2<dtype> {
    using type = tags::fusion;
    constexpr static std::size_t N_PRODUCTS = 1;
};

template<typename dtype>
struct Catalysis : public ReactionE2<dtype> {
    using type = tags::catalytic;
    constexpr static std::size_t N_PRODUCTS = 2;
};

}
}
