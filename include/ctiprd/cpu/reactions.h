//
// Created by mho on 4/13/22.
//

#pragma once

#include <tuple>
#include <vector>
#include <unordered_map>
#include <list>
#include <memory>

#include <spdlog/spdlog.h>
#include <tsl/robin_map.h>

#include <ctiprd/util/ops.h>
#include <ctiprd/util/hash.h>
#include <ctiprd/util/pbc.h>
#include <ctiprd/reactions/doi.h>

namespace ctiprd::cpu::reactions::impl {

template<typename Updater>
struct ReactionO1 {
    using dtype = typename Updater::dtype;
    using State = typename Updater::Position;
    [[nodiscard]] virtual bool shouldPerform([[maybe_unused]] const dtype &tau) const {
        return false;
    };

    virtual void operator()(std::size_t particleId, typename Updater::Particles &collection, Updater &updater) const {}
};

template<typename Updater>
struct ReactionO2 {
    using dtype = typename Updater::dtype;
    using State = typename Updater::Position;
    [[nodiscard]] virtual bool shouldPerform([[maybe_unused]] const dtype &tau) const { return false; };

    virtual void operator()(std::size_t id1, std::size_t id2, typename Updater::Particles &collection, Updater &updater) const {}

    dtype radiusSquared {};
    dtype rate {};
};

namespace detail {

template<typename dtype>
[[nodiscard]] bool shouldPerform(const dtype &tau, const dtype &rate) {
    return rnd::uniform_real<dtype>() < 1 - std::exp(-rate * tau);
}

template<typename Updater, typename Reaction>
struct CPUReactionO1 : ReactionO1<Updater> {};

template<typename Updater>
struct CPUReactionO1<Updater, ctiprd::reactions::doi::Decay<typename Updater::dtype>> : ReactionO1<Updater> {
    using ReactionType = ctiprd::reactions::doi::Decay<typename Updater::dtype>;
    using typename ReactionO1<Updater>::dtype;
    using typename ReactionO1<Updater>::State;

    explicit CPUReactionO1(const ReactionType &reaction) : baseReaction(reaction) {}

    bool shouldPerform(const dtype &tau) const override {
        return detail::shouldPerform(tau, baseReaction.rate);
    }

    void operator()(std::size_t id, typename Updater::Particles &collection, Updater &updater) const override {
        updater.remove(id, collection);
    }

    ReactionType baseReaction;
};

template<typename Updater>
struct CPUReactionO1<Updater, ctiprd::reactions::doi::Conversion<typename Updater::dtype>> : ReactionO1<Updater> {
    using ReactionType = ctiprd::reactions::doi::Conversion<typename Updater::dtype>;
    using typename ReactionO1<Updater>::dtype;
    using typename ReactionO1<Updater>::State;

    explicit CPUReactionO1(const ReactionType &reaction) : baseReaction(reaction) {}

    bool shouldPerform(const dtype &tau) const override {
        return detail::shouldPerform(tau, baseReaction.rate);
    }

    void operator()(std::size_t id, typename Updater::Particles &collection, Updater &updater) const override {
        updater.directUpdate(id, baseReaction.productType, {}, collection);
    }

    ReactionType baseReaction;
};

template<typename Updater>
struct CPUReactionO1<Updater, ctiprd::reactions::doi::Fission<typename Updater::dtype>> : ReactionO1<Updater> {
    using ReactionType = ctiprd::reactions::doi::Fission<typename Updater::dtype>;
    using typename ReactionO1<Updater>::dtype;
    using typename ReactionO1<Updater>::State;

    explicit CPUReactionO1(const ReactionType &reaction) : baseReaction(reaction) {}

    bool shouldPerform(const dtype &tau) const override {
        return detail::shouldPerform(tau, baseReaction.rate);
    }

    template<typename Generator>
    static auto normal(Generator &generator) {
        static std::normal_distribution<dtype> dist{0, 1};
        return dist(generator);
    }
    template<typename Generator>
    static auto uniform(Generator &generator) {
        static std::uniform_real_distribution<dtype> dist{0, 1};
        return dist(generator);
    }

    void operator()(std::size_t id, typename Updater::Particles &collection, Updater &updater) const override {
        const auto &c = collection.positionOf(id);
        State n {};
        std::transform(begin(n.data), end(n.data), begin(n.data), [](const auto &) {
            return normal(rnd::staticThreadLocalGenerator());
        });
        n /= n.norm();

        const auto distance = baseReaction.productDistance * std::pow(uniform(rnd::staticThreadLocalGenerator()), 1./Updater::dim);
        updater.add(baseReaction.productType2, c - 0.5 * distance * n, collection);
        updater.directUpdate(id, baseReaction.productType1, c + 0.5 * distance * n, collection);
    }

    ReactionType baseReaction;
};

template<typename Updater, typename Reaction>
struct CPUReactionO2 : ReactionO2<Updater> {};

template<typename Updater>
struct CPUReactionO2<Updater, ctiprd::reactions::doi::Fusion<typename Updater::dtype>> : public ReactionO2<Updater> {
    using ReactionType = ctiprd::reactions::doi::Fusion<typename Updater::dtype>;
    using Super = CPUReactionO2<Updater, ReactionType>;
    using typename ReactionO2<Updater>::dtype;
    using typename ReactionO2<Updater>::State;

    explicit CPUReactionO2(const ReactionType &reaction) : baseReaction(reaction) {
        Super::radiusSquared = baseReaction.reactionRadius * baseReaction.reactionRadius;
    }

    bool shouldPerform(const dtype &tau) const override {
        return detail::shouldPerform(tau, baseReaction.rate);
    }

    [[nodiscard]] auto key() const {
        return std::make_tuple(baseReaction.eductType1, baseReaction.eductType2);
    };

    void operator()(std::size_t id1, std::size_t id2, typename Updater::Particles &collection, Updater &updater) const override {
        if (updater.claim(id1)) {
            if (updater.claim(id2)) {
                auto w1 = collection.typeOf(id1) == baseReaction.eductType1 ? baseReaction.w1 : baseReaction.w2;
                updater.template directUpdate<false>(
                        id1, baseReaction.productType,
                        collection.positionOf(id1) + w1 * (collection.positionOf(id2) - collection.positionOf(id1)),
                        collection
                );
                updater.template remove<false>(id2, collection);
            } else {
                updater.release(id1);
            }
        }
    }

    ReactionType baseReaction;
};

template<typename Updater>
struct CPUReactionO2<Updater, ctiprd::reactions::doi::Catalysis<typename Updater::dtype>> : ReactionO2<Updater> {
    using ReactionType = ctiprd::reactions::doi::Catalysis<typename Updater::dtype>;
    using Super = CPUReactionO2<Updater, ReactionType>;
    using typename ReactionO2<Updater>::dtype;
    using typename ReactionO2<Updater>::State;

    explicit CPUReactionO2(const ReactionType &reaction) : baseReaction(reaction) {
        Super::radiusSquared = baseReaction.reactionRadius * baseReaction.reactionRadius;
    }

    bool shouldPerform(const dtype &tau) const override {
        return detail::shouldPerform(tau, baseReaction.rate);
    }

    [[nodiscard]] auto key() const {
        return std::make_tuple(baseReaction.catalyst, baseReaction.eductType);
    };

    void operator()(std::size_t id1, std::size_t id2, typename Updater::Particles &collection, Updater &updater) const override {
        if(collection.typeOf(id1) == baseReaction.catalyst) {
            updater.directUpdate(id2, baseReaction.productType, {}, collection);
        } else {
            updater.directUpdate(id1, baseReaction.productType, {}, collection);
        }
    }

    ReactionType baseReaction;
};

using KeyO2 = std::tuple<std::size_t, std::size_t>;
using HasherO2 = ctiprd::util::hash::ForwardBackwardTupleHasher<KeyO2>;
using EqO2 = ctiprd::util::hash::ForwardBackwardTupleEquality<KeyO2>;


}

template<typename Updater>
using ReactionsO1Map = tsl::robin_map<std::size_t, std::vector<const ReactionO1<Updater>*>>;
template<typename Updater>
using ReactionsO1Backing = std::list<std::unique_ptr<ReactionO1<Updater>>>;

template<typename Updater>
using ReactionsO2Map = tsl::robin_map<detail::KeyO2, std::vector<const ReactionO2<Updater>*>, detail::HasherO2, detail::EqO2>;
template<typename Updater>
using ReactionsO2Backing = std::list<std::unique_ptr<ReactionO2<Updater>>>;

template<typename Updater, typename System>
std::tuple<ReactionsO1Map<Updater>, ReactionsO1Backing<Updater>> generateMapO1(const System &system) {
    ReactionsO1Backing<Updater> backingData {};
    ReactionsO1Map<Updater> map;
    map.reserve(System::types.size());

    {
        for (const auto &[name, _] : System::types) {
            auto id = ctiprd::systems::particleTypeId<System::types>(name);
            map[id] = {};
        }
    }

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
std::tuple<ReactionsO2Map<Updater>, ReactionsO2Backing<Updater>> generateMapO2(const System &system) {
    ReactionsO2Backing<Updater> backingData {};
    ReactionsO2Map<Updater> map;
    map.reserve(ctiprd::util::binomialCoefficient(System::types.size(), 2));

    {
        for (const auto &[name1, d1] : System::types) {
            auto id1 = ctiprd::systems::particleTypeId<System::types>(name1);
            for (const auto &[name2, d2] : System::types) {
                auto id2 = ctiprd::systems::particleTypeId<System::types>(name2);
                map[{id1, id2}] = {};
            }
        }
    }

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
