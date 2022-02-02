//
// Created by mho on 5/5/21.
//

#pragma once

#include <algorithm>
#include <atomic>
#include <thread>

#include <ctiprd/util/ops.h>
#include <ctiprd/util/Index.h>
#include <ctiprd/thread/utils.h>
#include <ctiprd/ParticleCollection.h>

namespace ctiprd::nl {

namespace detail {
template<int dim, typename ParticleCollection, typename F>
void forEachParticle(ParticleCollection &collection, F op) {

    auto n = collection.nParticles();
    auto *forces = collection.forces().data();
    auto *pos = collection.positions();

    #pragma omp parallel for default(none), firstprivate(forces, pos, n, op)
    for (std::size_t i = 0; i < n; ++i) {
        op(i, pos + i * dim, forces + i * dim);
    }
}
}

template<int DIM, bool periodic, typename dtype>
class NeighborList {
public:
    NeighborList() : _gridSize() {}

    NeighborList(std::array<dtype, DIM> gridSize, dtype interactionRadius, std::size_t nParticles)
            : _gridSize(gridSize) {
        for (int i = 0; i < DIM; ++i) {
            _cellSize[i] = interactionRadius;
            if (gridSize[i] <= 0) throw std::invalid_argument("grid sizes must be positive.");
            nCells[i] = gridSize[i] / _cellSize[i];
        }
        _index = util::Index<DIM>(nCells);
        head.resize(_index.size());
        list.resize(nParticles + 1);
    }

    void update(ParticleCollection<DIM, dtype> &collection, int nJobs = 0) {
        std::fill(std::begin(list), std::end(list), 0);
        std::fill(std::begin(head), std::end(head), thread::copyable_atomic<std::size_t>());
        auto updateOp = [this](std::size_t particleId, dtype *pos, dtype * /*ignore*/) {
            auto boxId = positionToBoxIx(pos);

            // CAS
            auto &atomic = *head.at(boxId);
            auto currentHead = atomic.load();
            while (!atomic.compare_exchange_weak(currentHead, particleId)) {}
            list[particleId] = currentHead;
        };

        if (nJobs == 1) {
            for (std::size_t i = 0; i < collection.nParticles(); ++i) {
                updateOp(i, collection.positions() + DIM * i, nullptr);
            }
        } else {
            if (nJobs == 0) {
                nJobs = static_cast<decltype(nJobs)>(std::thread::hardware_concurrency());
            }
            if (nJobs <= 1) {
                throw std::logic_error("At this point nJobs should be >= 2");
            }
            std::size_t grainSize = collection.nParticles() / nJobs;
            auto *pptr = collection.positions();
            std::vector<std::jthread> jobs;
            for (int i = 0; i < nJobs - 1; ++i) {
                auto *pptrNext = std::min(pptr + grainSize * DIM,
                                          collection.positions() + collection.nParticles() * DIM);
                if (pptr != pptrNext) {
                    std::size_t idStart = i * grainSize;
                    jobs.emplace_back([&updateOp, idStart, pptr, pptrNext]() {
                        auto id = idStart;
                        for (auto *p = pptr; p != pptrNext; p += DIM, ++id) {
                            updateOp(id, p, nullptr);
                        }
                    });
                }
                pptr = pptrNext;
            }
            if (pptr != collection.positions() + DIM * collection.nParticles()) {
                auto pptrNext = collection.positions() + DIM * collection.nParticles();
                std::size_t idStart = (nJobs - 1) * grainSize;
                jobs.emplace_back([&updateOp, idStart, pptr, pptrNext]() {
                    auto id = idStart;
                    for (auto *p = pptr; p != pptrNext; p += DIM, ++id) {
                        updateOp(id, p, nullptr);
                    }
                });
            }
        }
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

    template<typename F>
    void forEachNeighbor(std::size_t id, ParticleCollection<DIM, dtype> &collection, F fun) const {
        auto *pos = collection.position(id);
        auto gridPos = this->gridPos(pos);
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
            auto boxId = positionToBoxIx(pos);
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
};

}
