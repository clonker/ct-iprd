//
// Created by mho on 4/13/22.
//

#pragma once

#include <tuple>
#include <vector>
#include <unordered_map>
#include <list>
#include <memory>

#include <ctiprd/util/hash.h>
#include <ctiprd/reactions/doi.h>
#include <spdlog/spdlog.h>

namespace ctiprd::cpu::reactions::doi {

template<typename Updater>
struct ReactionO1 {
    using dtype = typename Updater::dtype;
    using State = typename Updater::Position;
    [[nodiscard]] virtual bool shouldPerform(dtype tau, const State &s1, const std::size_t &t1) const { return false; };
};

template<typename Updater>
struct ReactionO2 {
    using dtype = typename Updater::dtype;
    using State = typename Updater::Position;
    [[nodiscard]] virtual bool shouldPerform(dtype tau, const State &s1, const std::size_t &t1,
                                             const State &s2, const std::size_t &t2) const { return false; };
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

template<typename Updater, typename Reaction>
struct CPUReactionO1 : ReactionO1<Updater> {};

template<typename Updater>
struct CPUReactionO1<Updater, ctiprd::reactions::doi::Decay<typename Updater::dtype>> : ReactionO1<Updater> {
    using ReactionType = ctiprd::reactions::doi::Decay<typename Updater::dtype>;
    using typename ReactionO1<Updater>::dtype;
    using typename ReactionO1<Updater>::State;

    explicit CPUReactionO1(const ReactionType &reaction) : baseReaction(reaction) {}

    bool shouldPerform(dtype tau, const State &s1, const std::size_t &t1) const override {
        return detail::shouldPerformO1(tau, s1, t1, baseReaction);
    }

    ReactionType baseReaction;
};

template<typename Updater>
struct CPUReactionO1<Updater, ctiprd::reactions::doi::Conversion<typename Updater::dtype>> : ReactionO1<Updater> {
    using ReactionType = ctiprd::reactions::doi::Conversion<typename Updater::dtype>;
    using typename ReactionO1<Updater>::dtype;
    using typename ReactionO1<Updater>::State;

    explicit CPUReactionO1(const ReactionType &reaction) : baseReaction(reaction) {}

    bool shouldPerform(dtype tau, const State &s1, const std::size_t &t1) const override {
        return detail::shouldPerformO1(tau, s1, t1, baseReaction);
    }

    ReactionType baseReaction;
};

template<typename Updater>
struct CPUReactionO1<Updater, ctiprd::reactions::doi::Fission<typename Updater::dtype>> : ReactionO1<Updater> {
    using ReactionType = ctiprd::reactions::doi::Fission<typename Updater::dtype>;
    using typename ReactionO1<Updater>::dtype;
    using typename ReactionO1<Updater>::State;

    explicit CPUReactionO1(const ReactionType &reaction) : baseReaction(reaction) {}

    bool shouldPerform(dtype tau, const State &s1, const std::size_t &t1) const override {
        return detail::shouldPerformO1(tau, s1, t1, baseReaction);
    }

    ReactionType baseReaction;
};

template<typename Updater, typename Reaction>
struct CPUReactionO2 : ReactionO2<Updater> {};

template<typename Updater>
struct CPUReactionO2<Updater, ctiprd::reactions::doi::Fusion<typename Updater::dtype>> : ReactionO2<Updater> {
    using ReactionType = ctiprd::reactions::doi::Fusion<typename Updater::dtype>;
    using typename ReactionO2<Updater>::dtype;
    using typename ReactionO2<Updater>::State;

    explicit CPUReactionO2(const ReactionType &reaction) : baseReaction(reaction) {}

    bool shouldPerform(dtype tau, const State &s1, const std::size_t &t1,
                       const State &s2, const std::size_t &t2) const override {
        return detail::shouldPerformO2(tau, s1, t1, s2, t2, baseReaction.eductType1, baseReaction.eductType2,
                                       baseReaction.rate, baseReaction.reactionRadius);
    }

    [[nodiscard]] auto key() const {
        return std::make_tuple(baseReaction.eductType1, baseReaction.eductType2);
    };

    ReactionType baseReaction;
};

template<typename Updater>
struct CPUReactionO2<Updater, ctiprd::reactions::doi::Catalysis<typename Updater::dtype>> : ReactionO2<Updater> {
    using ReactionType = ctiprd::reactions::doi::Catalysis<typename Updater::dtype>;
    using typename ReactionO2<Updater>::dtype;
    using typename ReactionO2<Updater>::State;

    explicit CPUReactionO2(const ReactionType &reaction) : baseReaction(reaction) {}

    bool shouldPerform(dtype tau, const State &s1, const std::size_t &t1,
                       const State &s2, const std::size_t &t2) const override {
        return detail::shouldPerformO2(tau, s1, t1, s2, t2, baseReaction.catalyst, baseReaction.eductType,
                                       baseReaction.rate, baseReaction.reactionRadius);
    }

    [[nodiscard]] auto key() const {
        return std::make_tuple(baseReaction.catalyst, baseReaction.eductType);
    };

    ReactionType baseReaction;
};

using KeyO2 = std::tuple<std::size_t, std::size_t>;
using HasherO2 = ctiprd::util::hash::ForwardBackwardTupleHasher<KeyO2>;
using EqO2 = ctiprd::util::hash::ForwardBackwardTupleEquality<KeyO2>;


}

template<typename Updater>
using ReactionsO1Map = std::unordered_map<std::size_t, std::vector<const ReactionO1<Updater>*>>;

template<typename Updater>
using ReactionsO2Map = std::unordered_map<detail::KeyO2, std::vector<const ReactionO2<Updater>*>, detail::HasherO2, detail::EqO2>;


template<typename Updater, typename System>
std::tuple<ReactionsO1Map<Updater>, std::list<std::unique_ptr<ReactionO1<Updater>>>> generateMapO1(const System &system) {
    std::list<std::unique_ptr<ReactionO1<Updater>>> backingData {};
    ReactionsO1Map<Updater> map;

    [&]<auto... I>(std::index_sequence<I...>) {
        ([&](const auto &reaction) {
            using ReactionType = std::tuple_element_t<I, typename System::ReactionsO1>;
            backingData.push_back(std::make_unique<detail::CPUReactionO1<Updater, ReactionType>>(reaction));
            const auto &ref = backingData.back();
            map[reaction.eductType].push_back(ref.get());
        }(std::get<I>(system.reactionsO1)), ...);
    }(std::make_index_sequence<std::tuple_size_v<typename System::ReactionsO1>>{});
    return std::make_tuple(map, std::move(backingData));
}

template<typename Updater, typename System>
std::tuple<ReactionsO2Map<Updater>, std::list<std::unique_ptr<ReactionO2<Updater>>>> generateMapO2(const System &system) {
    std::list<std::unique_ptr<ReactionO2<Updater>>> backingData {};
    ReactionsO2Map<Updater> map;

    [&]<auto... I>(std::index_sequence<I...>) {
        ([&](const auto &reaction) {
            using ReactionType = std::tuple_element_t<I, typename System::ReactionsO2>;
            using ReactionImpl = detail::CPUReactionO2<Updater, ReactionType>;
            backingData.push_back(std::make_unique<ReactionImpl>(reaction));
            const auto &ref = backingData.back();
            map[ReactionImpl{reaction}.key()].push_back(ref.get());
        }(std::get<I>(system.reactionsO2)), ...);
    }(std::make_index_sequence<std::tuple_size_v<typename System::ReactionsO2>>{});
    return std::make_tuple(map, std::move(backingData));
}

}
