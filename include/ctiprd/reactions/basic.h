/**
 *
 *
 * @file basic.h
 * @brief 
 * @author clonker
 * @date 3/27/22
 */
#pragma once

namespace ctiprd::reactions {

template<typename dtype, std::size_t eductId, std::size_t... productIds >
struct ReactionE1 {
    constexpr static std::size_t N_EDUCTS = 1;
    constexpr static std::size_t N_PRODUCTS = sizeof...(productIds);
    constexpr static std::size_t particleType = eductId;
    constexpr static std::array<std::size_t, N_PRODUCTS> productTypes {productIds...};

    template<typename State, typename ParticleType>
    [[nodiscard]] bool shouldPerform(dtype tau, const State &state, const ParticleType &t) const {
        return t == eductId && rnd::uniform_real<dtype>() < 1 - std::exp(-rate * tau);
    }

    dtype rate;
    dtype productRadius;
};

template<typename dtype, std::size_t eductId1, std::size_t eductId2, std::size_t... productIds >
struct ReactionE2 {
    constexpr static std::size_t N_EDUCTS = 2;
    constexpr static std::size_t N_PRODUCTS = sizeof...(productIds);
    constexpr static std::size_t particleType1 = eductId1;
    constexpr static std::size_t particleType2 = eductId2;
    constexpr static std::array<std::size_t, N_PRODUCTS> productTypes {productIds...};

    template<typename State, typename ParticleType>
    [[nodiscard]] bool shouldPerform(dtype tau, const State &state, const ParticleType &t,
                                     const State &state2, const ParticleType &t2) const {
        if ((t == eductId1 && t2 == eductId2) || (t == eductId2 && t2 == eductId1)) {
            return (state - state2).normSquared() < productRadius * productRadius && rnd::uniform_real<dtype>() < 1 - std::exp(-rate * tau);
        }
        return false;
    }

    dtype rate;
    dtype productRadius;
};

template<typename dtype, std::size_t id>
using Decay = ReactionE1<dtype, id>;

template<typename dtype, std::size_t id1, std::size_t id2>
using Conversion = ReactionE1<dtype, id1, id2>;

template<typename dtype, std::size_t id1, std::size_t id2, std::size_t id3>
using Fission = ReactionE1<dtype, id1, id2, id3>;

template<typename dtype, std::size_t id1, std::size_t id2, std::size_t id3>
using Fusion = ReactionE2<dtype, id1, id2, id3>;

template<typename dtype, std::size_t catalyst, std::size_t id1, std::size_t id2>
using Catalysis = ReactionE2<dtype, catalyst, id1, catalyst, id2>;
}
