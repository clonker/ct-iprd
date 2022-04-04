/**
 * @file UncontrolledApproximation.h
 * @brief 
 * @author clonker
 * @date 3/27/22
 */
#pragma once

#include <variant>
#include <execution>
#include <algorithm>

#include <ctiprd/any.h>
#include <ctiprd/NeighborList.h>
#include "basic.h"

namespace ctiprd::reactions {

namespace impl {

template<typename Reaction>
struct ReactionImpl;

template<typename dtype>
struct ReactionImpl<doi::Decay<dtype>> {
    explicit ReactionImpl(const doi::Decay<dtype> *reaction) {}

    template<typename Updater>
    void operator()(auto id1, auto /*ignore*/, typename Updater::Particles &collection, Updater &updater) {
        updater.remove(id1, collection);
    }
};

template<typename dtype>
struct ReactionImpl<doi::Conversion<dtype>> {
    explicit ReactionImpl(const doi::Conversion<dtype> *reaction) : reaction(reaction) {}

    template<typename Updater>
    void operator()(auto id1, auto /*ignore*/, typename Updater::Particles &collection, Updater &updater) {
        updater.directUpdate(id1, reaction->productType, {}, collection);
    }

    const doi::Conversion<dtype> *reaction;
};

template<typename dtype>
struct ReactionImpl<doi::Catalysis<dtype>> {
    explicit ReactionImpl(const doi::Catalysis<dtype> *reaction) : reaction(reaction) {}

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

}

template<typename dtype, typename Updater>
struct ReactionEvent {
    std::size_t id1, id2;
    std::function<void(std::size_t id1, std::size_t id2, typename Updater::Particles &, Updater&)> perform;

    void operator()(typename Updater::Particles &particles, Updater &updater) {
        perform(id1, id2, particles, updater);
    }
};

template<typename System, typename Generator>
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

        using Updater = ParticleCollectionUpdater<Particles>;

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

                std::apply([&pos, &type, &events, &tau, &id](auto &&... reaction) {
                    (
                            ([&events, &pos, &type, &tau, &id](const auto &r) {
                                if (r.shouldPerform(tau, pos, type)) {
                                    events.push_back(ReactionEvent<dtype, Updater>{
                                            .id1 = id,
                                            .id2 = id,
                                            .perform = impl::ReactionImpl<std::decay_t<decltype(r)>>{&r}
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
                                        if (r.shouldPerform(tau, pos, type, neighborPos, neighborType)) {
                                            events.push_back(ReactionEvent<dtype, Updater>{
                                                    .id1 = id,
                                                    .id2 = neighborId,
                                                    .perform = impl::ReactionImpl<std::decay_t<decltype(r)>>{&r}
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

        {
            std::shuffle(begin(events), end(events), rnd::staticThreadLocalGenerator<Generator>());
            ParticleCollectionUpdater updater {*particles};

            for(auto it = begin(events); it != end(events); ++it) {
                (*it)(*particles, updater);
            }

            particles->update(begin(updater.toAdd), end(updater.toAdd), begin(updater.toRemove), end(updater.toRemove));
        }
    }

    std::unique_ptr<NeighborList> neighborList_;
};

}
