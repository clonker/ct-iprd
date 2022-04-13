//
// Created by mho on 4/13/22.
//

#pragma once

#include <tuple>
#include <vector>
#include <unordered_map>

#include <ctiprd/util/hash.h>
#include <ctiprd/reactions/doi.h>

namespace ctiprd::cpu::reactions::doi {

template<typename ParticleCollection>
struct ReactionO1 {
    using dtype = typename ParticleCollection::dtype;
    using State = typename ParticleCollection::Position;

    [[nodiscard]] virtual bool shouldPerform(dtype tau, const State &s1, const std::size_t &t1) const = 0;
};

template<typename ParticleCollection>
struct ReactionO2 {
    using dtype = typename ParticleCollection::dtype;
    using State = typename ParticleCollection::Position;

    [[nodiscard]] virtual bool shouldPerform(dtype tau, const State &s1, const std::size_t &t1,
                                             const State &s2, const std::size_t &t2) const = 0;
};

namespace detail {

template<typename dtype, typename State, typename ParticleType, typename Reaction>
[[nodiscard]] bool shouldPerformO1(dtype tau, const State &state, const ParticleType &t, const Reaction &reaction) {
    return t == reaction.eductType && rnd::uniform_real<dtype>() < 1 - std::exp(-reaction.rate * tau);
}

template<typename dtype, typename State, typename ParticleType>
[[nodiscard]] bool shouldPerformO2(dtype tau, const State &state, const ParticleType &t,
                                   const State &state2, const ParticleType &t2,
                                   const ParticleType &eductType1, const ParticleType &eductType2,
                                   dtype rate, dtype radius) {
    if ((t == eductType1 && t2 == eductType2) || (t == eductType2 && t2 == eductType1)) {
        return (state - state2).normSquared() < radius * radius &&
               rnd::uniform_real<dtype>() < 1 - std::exp(-rate * tau);
    }
    return false;
}

template<typename ParticleCollection, typename Reaction>
struct CPUReaction;

template<typename ParticleCollection>
struct CPUReaction<ParticleCollection, ctiprd::reactions::doi::Decay<typename ParticleCollection::dtype>>
        : ReactionO1<ParticleCollection>, ctiprd::reactions::doi::Decay<typename ParticleCollection::dtype> {
    using SuperReaction = ctiprd::reactions::doi::Decay<typename ParticleCollection::dtype>;
    using typename ReactionO1<ParticleCollection>::dtype;
    using typename ReactionO1<ParticleCollection>::State;

    explicit CPUReaction(const SuperReaction &reaction) : SuperReaction(reaction) {}

    bool shouldPerform(dtype tau, const State &s1, const std::size_t &t1) const override {
        return detail::shouldPerformO1(tau, s1, t1, *this);
    }
};

using KeyO2 = std::tuple<std::size_t, std::size_t>;
using HasherO2 = ctiprd::util::hash::ForwardBackwardTupleHasher<KeyO2>;
using EqO2 = ctiprd::util::hash::ForwardBackwardTupleEquality<KeyO2>;


}

template<typename ParticleCollection>
using ReactionsO1Map = std::unordered_map<std::size_t, std::vector<ReactionO1<ParticleCollection>*>>;

template<typename ParticleCollection>
using ReactionsO2Map = std::unordered_map<detail::KeyO2, std::vector<ReactionO2<ParticleCollection>*>, detail::HasherO2, detail::EqO2>;


auto generateMapO1() {
}

}
