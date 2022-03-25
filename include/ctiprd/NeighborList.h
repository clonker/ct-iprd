//
// Created by mho on 5/5/21.
//

#pragma once

#include <algorithm>
#include <atomic>
#include <thread>

#include <ctiprd/config.h>
#include <ctiprd/util/ops.h>
#include <ctiprd/util/Index.h>
#include <ctiprd/thread/utils.h>

namespace ctiprd::nl {

template<int DIM, bool periodic, typename dtype, typename Pool=config::ThreadPool>
class NeighborList {
public:
    NeighborList(std::array<dtype, DIM> gridSize, dtype interactionRadius, std::shared_ptr<Pool> pool)
            : _gridSize(gridSize), pool(pool) {
        for (int i = 0; i < DIM; ++i) {
            _cellSize[i] = interactionRadius;
            if (gridSize[i] <= 0) throw std::invalid_argument("grid sizes must be positive.");
            nCells[i] = gridSize[i] / _cellSize[i];
        }
        _index = util::Index<DIM>(nCells);
        head.resize(_index.size());
    }

    ~NeighborList() = default;
    NeighborList(const NeighborList&) = delete;
    NeighborList &operator=(const NeighborList&) = delete;

    NeighborList(NeighborList&&) noexcept = default;
    NeighborList &operator=(NeighborList&&) noexcept = default;

    template<typename ParticleCollection>
    void update(std::shared_ptr<ParticleCollection> collection) {
        list.resize(collection->nParticles() + 1);
        std::fill(std::begin(list), std::end(list), 0);
        std::fill(std::begin(head), std::end(head), thread::copyable_atomic<std::size_t>());
        auto updateOp = [this](std::size_t particleId, const auto &pos, const auto &/*ignore*/) {
            auto boxId = positionToBoxIx(&pos.data[0]);

            // CAS
            auto &atomic = *head.at(boxId);
            auto currentHead = atomic.load();
            while (!atomic.compare_exchange_weak(currentHead, particleId)) {}
            list[particleId] = currentHead;
        };

        collection->template forEachParticle(updateOp, pool);
        pool->waitForTasks();
    }

    typename util::Index<DIM>::GridDims gridPos(const dtype *pos) const {
        std::array<std::uint32_t, DIM> projections;
        for (auto i = 0u; i < DIM; ++i) {
            projections[i] = static_cast<std::uint32_t>(
                    std::max(static_cast<dtype>(0.), static_cast<dtype>(std::floor((pos[i] + .5 * _gridSize[i]) / _cellSize[i]))));
            projections[i] = util::clamp(projections[i], 0U, nCells[i] - 1);
        }
        return projections;
    }

    std::uint32_t positionToBoxIx(const dtype *pos) const {
        auto boxId = _index.index(gridPos(pos));
        return boxId;
    }

    template<typename ParticleCollection, typename F>
    void forEachNeighbor(std::size_t id, ParticleCollection &collection, F &&fun) const {
        const auto &pos = collection.position(id);
        auto gridPos = this->gridPos(&pos[0]);
        for (int i = 0u; i < DIM; ++i) {
            std::uint32_t begin, end;
            if constexpr(!periodic) {
                begin = gridPos[i] > 0 ? (gridPos[i] - 1) : gridPos[i];
                end = gridPos[i] < nCells[i] - 1 ? (gridPos[i] + 2) : gridPos[i] + 1;
            } else {
                begin = gridPos[i] > 0 ? gridPos[i] - 1 : nCells[i] - 1;
                end = (gridPos[i] + 2) % nCells[i];
            }

            auto cellPos = gridPos;
            for (auto k = begin; k != end; k = periodic ? (k + 1) % nCells[i] : k + 1) {
                cellPos[i] = k;

                if (cellPos != gridPos) {
                    auto cellId = _index.index(cellPos);
                    auto neighborId = (*head.at(cellId)).load();
                    while (neighborId != 0) {
                        fun(neighborId, collection.position(neighborId), collection.force(neighborId));
                        neighborId = list.at(neighborId);
                    }
                }
            }
        }
        {
            auto boxId = positionToBoxIx(&pos[0]);
            auto neighborId = (*head.at(boxId)).load();
            while (neighborId != 0) {
                if (neighborId != id) {
                    fun(neighborId, collection.position(neighborId), collection.force(neighborId));
                }
                neighborId = list.at(neighborId);
            }
        }
    }

    const std::array<dtype, DIM> &gridSize() const {
        return _gridSize;
    }

private:
    std::array<dtype, DIM> _cellSize{};
    std::array<dtype, DIM> _gridSize{};
    std::vector<thread::copyable_atomic<std::size_t>> head{};
    std::vector<std::size_t> list{};
    util::Index<DIM> _index{};
    std::array<std::uint32_t, DIM> nCells{};

    std::shared_ptr<Pool> pool;
};

}
