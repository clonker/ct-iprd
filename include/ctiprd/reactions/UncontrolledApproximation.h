/**
 * @file UncontrolledApproximation.h
 * @brief 
 * @author clonker
 * @date 3/27/22
 */
#pragma once

#include <variant>
#include <algorithm>

#include <ctiprd/any.h>
#include <ctiprd/NeighborList.h>
#include <ctiprd/reactions/basic.h>
#include <ctiprd/ParticleCollection.h>

namespace ctiprd::reactions {

namespace impl {

namespace detail {

template<typename dtype, typename State, typename ParticleType, typename Reaction>
[[nodiscard]] bool shouldPerformO1(dtype tau, const State &state, const ParticleType &t, const Reaction* reaction) {
    return t == reaction->eductType && rnd::uniform_real<dtype>() < 1 - std::exp(-reaction->rate * tau);
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

}

template<typename Reaction>
struct ReactionImpl;

template<typename dtype>
struct ReactionImpl<doi::Decay<dtype>> {
    explicit ReactionImpl(const doi::Decay<dtype> *reaction) : reaction(reaction) {}

    template<typename State, typename ParticleType>
    [[nodiscard]] bool shouldPerform(dtype tau, const State &state, const ParticleType &t) const {
        return detail::shouldPerformO1(tau, state, t, reaction);
    }

    template<typename Updater>
    void operator()(auto id1, auto /*ignore*/, typename Updater::Particles &collection, Updater &updater) {
        updater.remove(id1, collection);
    }

    const doi::Decay<dtype> *reaction;
};

template<typename dtype>
struct ReactionImpl<doi::Conversion<dtype>> {
    explicit ReactionImpl(const doi::Conversion<dtype> *reaction) : reaction(reaction) {}

    template<typename State, typename ParticleType>
    [[nodiscard]] bool shouldPerform(dtype tau, const State &state, const ParticleType &t) const {
        return detail::shouldPerformO1(tau, state, t, reaction);
    }

    template<typename Updater>
    void operator()(auto id1, auto /*ignore*/, typename Updater::Particles &collection, Updater &updater) {
        updater.directUpdate(id1, reaction->productType, {}, collection);
    }

    const doi::Conversion<dtype> *reaction;
};

template<typename dtype>
struct ReactionImpl<doi::Fission<dtype>> {
    explicit ReactionImpl(const doi::Fission<dtype> *reaction) : reaction(reaction) {}

    template<typename State, typename ParticleType>
    [[nodiscard]] bool shouldPerform(dtype tau, const State &state, const ParticleType &t) const {
        return detail::shouldPerformO1(tau, state, t, reaction);
    }

    template<typename Updater, typename Particles = typename Updater::Particles>
    void operator()(auto id1, auto /*ignore*/, Particles &collection, Updater &updater) {
        const auto &c = collection.positionOf(id1);
        typename Particles::Position n {};
        std::transform(begin(n.data), end(n.data), begin(n.data), [](const auto &) {
            return normal(rnd::staticThreadLocalGenerator());
        });
        n /= n.norm();

        const auto distance = reaction->productDistance * std::pow(uniform(rnd::staticThreadLocalGenerator()), 1/Particles::DIM);
        updater.add(reaction->productType2, c - 0.5 * distance * n);
        updater.directUpdate(id1, reaction->productType1, c + 0.5 * distance * n, collection);
    }

    const doi::Fission<dtype> *reaction;

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
};

template<typename dtype>
struct ReactionImpl<doi::Catalysis<dtype>> {
    explicit ReactionImpl(const doi::Catalysis<dtype> *reaction) : reaction(reaction) {}

    template<typename State, typename ParticleType>
    [[nodiscard]] bool shouldPerform(dtype tau, const State &state, const ParticleType &t,
                                     const State &state2, const ParticleType &t2) const {
        return detail::shouldPerformO2(tau, state, t, state2, t2,
                                       reaction->catalyst, reaction->eductType, reaction->rate,
                                       reaction->reactionRadius);
    }

    template<typename Updater>
    void operator()(auto id1, auto id2, typename Updater::Particles &collection, Updater &updater) {
        if(collection.typeOf(id1) == reaction->catalyst) {
            updater.directUpdate(id2, reaction->productType, {}, collection);
        } else {
            updater.directUpdate(id1, reaction->productType, {}, collection);
        }
    }

    const doi::Catalysis<dtype> *reaction;
};

template<typename dtype>
struct ReactionImpl<doi::Fusion<dtype>> {
    explicit ReactionImpl(const doi::Fusion<dtype> *reaction) : reaction(reaction) {}

    template<typename State, typename ParticleType>
    [[nodiscard]] bool shouldPerform(dtype tau, const State &state, const ParticleType &t,
                                     const State &state2, const ParticleType &t2) const {
        return detail::shouldPerformO2(tau, state, t, state2, t2, reaction->eductType1, reaction->eductType2,
                                       reaction->rate, reaction->reactionRadius);
    }

    template<typename Updater>
    void operator()(auto id1, auto id2, typename Updater::Particles &collection, Updater &updater) {
        if (updater.claim(id1)) {
            if (updater.claim(id2)) {
                auto w1 = collection.typeOf(id1) == reaction->eductType1 ? reaction->w1 : reaction->w2;
                updater.template directUpdate<false>(
                        id1, reaction->productType,
                        collection.positionOf(id1) + w1 * (collection.positionOf(id2) - collection.positionOf(id1)),
                        collection
                );
                updater.template remove<false>(id2, collection);
            } else {
                updater.release(id1);
            }
        }
    }
    const doi::Fusion<dtype> *reaction;
};

}

template<typename dtype, typename Updater>
struct ReactionEvent {
    std::size_t id1, id2;
    std::function<void(std::size_t id1, std::size_t id2, typename Updater::Particles &, Updater&)> perform;

    void operator()(typename Updater::Particles &particles, Updater &updater) {
        perform(id1, id2, particles, updater);
    }
};

template<typename System, typename Generator=config::DefaultGenerator>
struct UncontrolledApproximation {
    static constexpr std::size_t DIM = System::DIM;
    using dtype = typename System::dtype;
    using NeighborList = ctiprd::nl::NeighborList<DIM, System::periodic, dtype>;

    using ReactionsO1 = typename System::ReactionsO1;
    using ReactionsO2 = typename System::ReactionsO2;


    static constexpr int nReactionsO1 = std::tuple_size_v<ReactionsO1>;
    static constexpr int nReactionsO2 = std::tuple_size_v<ReactionsO2>;

    template<typename Particles, typename Pool>
    void reactions(const System &system, dtype tau, std::shared_ptr<Particles> particles,
                   std::shared_ptr<Pool> pool) {

        using Updater = ParticleCollectionUpdater<System, Particles>;

        if constexpr(nReactionsO2 > 0) {
            if (!neighborList_) {
                auto cutoff = reactions::reactionRadius<dtype>(system.reactionsO2);
                neighborList_ = std::make_unique<NeighborList>(System::boxSize, cutoff);
            }
            neighborList_->update(particles, pool);
        }

        std::mutex mutex;
        std::vector<ReactionEvent<dtype, Updater>> events;
        {
            auto worker = [&system, tau, nl = neighborList_.get(), data = particles.get(), &events, &mutex](
                    auto id, typename Particles::Position &pos, const typename Particles::ParticleType &type,
                    typename Particles::Force &force
            ) {
                std::vector<ReactionEvent<dtype, Updater>> localEvents;

                std::apply([&pos, &type, &localEvents, &tau, &id](auto &&... reaction) {
                    (
                            ([&localEvents, &pos, &type, &tau, &id](const auto &r) {
                                impl::ReactionImpl<std::decay_t<decltype(r)>> impl {&r};
                                if (impl.shouldPerform(tau, pos, type)) {
                                    localEvents.push_back(ReactionEvent<dtype, Updater>{
                                            .id1 = id,
                                            .id2 = id,
                                            .perform = impl
                                    });
                                }
                            }(reaction)),
                            ...);
                }, system.reactionsO1);

                if constexpr(nReactionsO2 > 0) {
                    nl->template forEachNeighbor(id, *data, [&](auto neighborId, const auto &neighborPos,
                                                                const auto &neighborType,
                                                                const auto &neighborForce) {
                        std::apply([&](auto &&... reaction) {
                            (
                                    ([&](const auto &r) {
                                        impl::ReactionImpl<std::decay_t<decltype(r)>> impl {&r};
                                        if (impl.shouldPerform(tau, pos, type, neighborPos, neighborType)) {
                                            localEvents.push_back(ReactionEvent<dtype, Updater>{
                                                    .id1 = id,
                                                    .id2 = neighborId,
                                                    .perform = impl
                                            });
                                        }
                                    }(reaction)),
                                    ...);
                        }, system.reactionsO2);
                    });
                }

                {
                    std::scoped_lock lock{mutex};
                    events.reserve(events.size() + localEvents.size());
                    events.insert(end(events), begin(localEvents), end(localEvents));
                }
            };

            particles->forEachParticle(worker, pool);
        }
        pool->waitForTasks();
        spdlog::debug("got reactions {}", events.size());

        {
            std::shuffle(begin(events), end(events), rnd::staticThreadLocalGenerator<Generator>());
            Updater updater {*particles};

            for(auto it = begin(events); it != end(events); ++it) {
                (*it)(*particles, updater);
            }

            particles->update(begin(updater.toAdd), end(updater.toAdd), begin(updater.toRemove), end(updater.toRemove));
        }
    }

    std::unique_ptr<NeighborList> neighborList_;
};

}
