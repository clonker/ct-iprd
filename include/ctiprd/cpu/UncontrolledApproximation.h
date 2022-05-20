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
    std::uint8_t nEducts {};
    std::size_t id1 {}, id2 {};
    std::size_t reactionIndex {};
    std::size_t type1 {}, type2 {};
    bool valid {true};
};

template<typename ParticleCollection, typename System, typename Generator=config::DefaultGenerator>
struct UncontrolledApproximation {
    static constexpr std::size_t DIM = System::DIM;
    using dtype = typename System::dtype;
    using NeighborList = nl::NeighborList<DIM, System::periodic, dtype, false>;

    using ReactionsO1 = typename System::ReactionsO1;
    using ReactionsO2 = typename System::ReactionsO2;

    using Updater = SynchronizedParticleCollectionUpdater<System, ParticleCollection>;

    static constexpr int nReactionsO1 = std::tuple_size_v<ReactionsO1>;
    static constexpr int nReactionsO2 = std::tuple_size_v<ReactionsO2>;

    explicit UncontrolledApproximation(const System &system) {
        if constexpr(nReactionsO2 > 0) {
            auto cutoff = ctiprd::reactions::reactionRadius<dtype>(system.reactionsO2);
            neighborList_ = std::make_unique<NeighborList>(System::boxSize, cutoff);
        }

        std::tie(reactionsO1, backingO1) = reactions::impl::generateMapO1<Updater>(system);
        std::tie(reactionsO2, backingO2) = reactions::impl::generateMapO2<Updater>(system);

        if(neighborList_) {
            std::unordered_set<std::size_t> neighborListTypes{};
            for (const auto &reactionElement: reactionsO2) {
                if (!reactionElement.second.empty()) {
                    neighborListTypes.emplace(std::get<0>(reactionElement.first));
                    neighborListTypes.emplace(std::get<1>(reactionElement.first));
                }
            }
            neighborList_->setTypes(neighborListTypes);
        }
    }

    template<typename Pool>
    void reactions(const dtype tau, std::shared_ptr<ParticleCollection> particles, std::shared_ptr<Pool> pool) {

        std::vector<std::future<void>> futures;

        if constexpr(nReactionsO2 > 0) {
            neighborList_->update(particles, pool);
        }

        std::mutex mutex;
        std::vector<ReactionEvent> events;
        {
            {
                const auto worker = [this, tau, neighborList = neighborList_.get(), data = particles.get(), &events, &mutex](
                        const auto particleId, typename ParticleCollection::Position &pos,
                        const typename ParticleCollection::ParticleType &type,
                        const auto &/*ignore*/
                ) {
                    std::vector<ReactionEvent> localEvents;
                    const auto &reactions = reactionsO1[type];
                    for (std::size_t i = 0; i < reactions.size(); ++i) {
                        if (reactions[i]->shouldPerform(tau)) {
                            localEvents.push_back({1, particleId, particleId, i, type, 0});
                        }
                    }

                    {
                        const std::scoped_lock lock{mutex};
                        events.reserve(events.size() + localEvents.size());
                        events.insert(end(events), begin(localEvents), end(localEvents));
                    }
                };

                auto fun = particles->forEachParticle(worker, pool);
                std::move(begin(fun), end(fun), std::back_inserter(futures));
            }

            if constexpr(nReactionsO2 > 0) {
                const auto that = this;
                const auto worker = [that, tau, &data = *particles, &mutex, &events](const auto &cellIndex) {
                    std::vector<ReactionEvent> localEvents;

                    that->neighborList_->template forEachNeighborInCell<false>([that, tau, &localEvents, &data](const auto &id1, const auto &id2) {
                        const auto type1 = data.typeOf(id1);
                        const auto type2 = data.typeOf(id2);
                        const auto &position1 = data.positionOf(id1);
                        const auto &position2 = data.positionOf(id2);

                        const auto &reactions = that->reactionsO2[{type1, type2}];

                        const auto distsq = util::pbc::shortestDifference<System>(position1, position2).normSquared();
                        for (std::size_t i = 0; i < reactions.size(); ++i) {
                            if (distsq <= reactions[i]->radiusSquared && reactions[i]->shouldPerform(tau)) {
                                localEvents.push_back({2, id1, id2, i, type1, type2});
                            }
                        }
                    }, cellIndex);

                    {
                        std::scoped_lock lock{mutex};
                        events.reserve(events.size() + localEvents.size());
                        events.insert(end(events), begin(localEvents), end(localEvents));
                    }
                };
                auto cellFutures = neighborList_->forEachCell(worker, pool);
                std::move(begin(cellFutures), end(cellFutures), std::back_inserter(futures));
            }
        }
        for (auto &future : futures) { future.wait(); }

        {
            std::shuffle(begin(events), end(events), rnd::staticThreadLocalGenerator<Generator>());
            Updater updater {*particles};

            for(auto it = begin(events); it != end(events); ++it) {
                if(it->valid) {
                    if (it->nEducts == 1) {
                        (*reactionsO1[it->type1][it->reactionIndex])(it->id1, *particles, updater);
                    } else {
                        (*reactionsO2[{it->type1, it->type2}][it->reactionIndex])(it->id1, it->id2, *particles, updater);
                    }
                    for (auto it2 = it + 1; it2 != end(events); ++it2) {
                        if(it2->valid && (it->id1 == it2->id1 || it->id1 == it2->id2 || it->id2 == it2->id1 || it->id2 == it2->id2)) {
                            it2->valid = false;
                        }
                    }
                }
            }
        }
        if constexpr(nReactionsO2 > 0) {
            // particles->sort();
        }
    }

    std::unique_ptr<NeighborList> neighborList_;
    reactions::impl::ReactionsO1Map<Updater> reactionsO1;
    reactions::impl::ReactionsO1Backing<Updater> backingO1;
    reactions::impl::ReactionsO2Map<Updater> reactionsO2;
    reactions::impl::ReactionsO2Backing<Updater> backingO2;
};

}
