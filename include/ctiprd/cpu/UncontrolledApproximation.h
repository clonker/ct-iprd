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
#include <ctiprd/reactions/doi.h>
#include <ctiprd/cpu/NeighborList.h>
#include <ctiprd/cpu/ParticleCollection.h>
#include <ctiprd/cpu/reactions.h>

namespace ctiprd::cpu {

struct ReactionEvent {
    std::uint8_t nEducts;
    std::size_t id1, id2;
    std::size_t reactionIndex;
};

template<typename ParticleCollection, typename System, typename Generator=config::DefaultGenerator>
struct UncontrolledApproximation {
    static constexpr std::size_t DIM = System::DIM;
    using dtype = typename System::dtype;
    using NeighborList = nl::NeighborList<DIM, System::periodic, dtype>;

    using ReactionsO1 = typename System::ReactionsO1;
    using ReactionsO2 = typename System::ReactionsO2;

    using Updater = ParticleCollectionUpdater<System, ParticleCollection>;

    static constexpr int nReactionsO1 = std::tuple_size_v<ReactionsO1>;
    static constexpr int nReactionsO2 = std::tuple_size_v<ReactionsO2>;

    explicit UncontrolledApproximation(const System &system) {
        if constexpr(nReactionsO2 > 0) {
            if (!neighborList_) {
                auto cutoff = ctiprd::reactions::reactionRadius<dtype>(system.reactionsO2);
                neighborList_ = std::make_unique<NeighborList>(System::boxSize, cutoff);
            }
        }

        std::tie(reactionsO1, backingO1) = reactions::impl::generateMapO1<Updater>(system);
        std::tie(reactionsO2, backingO2) = reactions::impl::generateMapO2<Updater>(system);
    }

    template<typename Pool>
    void reactions(const System &system, dtype tau, std::shared_ptr<ParticleCollection> particles,
                   std::shared_ptr<Pool> pool) {

        if constexpr(nReactionsO2 > 0) {
            neighborList_->update(particles, pool);
        }

        std::mutex mutex;
        std::vector<ReactionEvent> events;
        {
            auto worker = [&system, tau, nl = neighborList_.get(), data = particles.get(), &events, &mutex](
                    auto id, typename ParticleCollection::Position &pos, const typename ParticleCollection::ParticleType &type,
                    typename ParticleCollection::Force &force
            ) {
                std::vector<ReactionEvent> localEvents;

                /*std::apply([&pos, &type, &localEvents, &tau, &id](auto &&... reaction) {
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
                }*/

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
            Updater updater {*particles};

            for(auto it = begin(events); it != end(events); ++it) {
                //(*it)(*particles, updater);
            }

            particles->update(begin(updater.toAdd), end(updater.toAdd), begin(updater.toRemove), end(updater.toRemove));
        }
    }

    std::unique_ptr<NeighborList> neighborList_;
    reactions::impl::ReactionsO1Map<Updater> reactionsO1;
    reactions::impl::ReactionsO1Backing<Updater> backingO1;
    reactions::impl::ReactionsO2Map<Updater> reactionsO2;
    reactions::impl::ReactionsO2Backing<Updater> backingO2;
};

}
