//
// Created by mho on 5/5/21.
//

#pragma once

#include <algorithm>
#include <atomic>
#include <thread>
#include <unordered_set>

#include <ctiprd/config.h>
#include <ctiprd/util/ops.h>
#include <ctiprd/util/Index.h>
#include <ctiprd/thread/utils.h>

namespace ctiprd::cpu::nl {

namespace detail {

template<bool periodic, typename Index>
struct CellAdjacency {
    CellAdjacency(const Index &index, std::size_t r) {
        typename Index::GridDims nNeighbors {};
        for (std::size_t i = 0; i < Index::Dims; ++i) {
            nNeighbors[i] = std::min(index[i], static_cast<std::size_t>(2 * r + 1));
        }
        auto nAdjacentCells = std::accumulate(begin(nNeighbors), end(nNeighbors), 1, std::multiplies<>());
        cellNeighbors = util::Index<2>{std::array<std::size_t, 2>{index.size(), 1 + nAdjacentCells}};
        cellNeighborsContent.resize(cellNeighbors.size());

        for(std::size_t i = 0; i < index.size(); ++i) {
            auto ijk = index.inverse(i);

        }
    }

    // index of size (n_cells x (1 + nAdjacentCells)), where the first element tells how many adj cells are stored
    util::Index<2> cellNeighbors;
    // backing vector of _cellNeighbors index of size (n_cells x (1 + nAdjacentCells))
    std::vector<std::size_t> cellNeighborsContent;
};
}

template<int DIM, bool periodic, typename dtype, bool allTypes=true>
class NeighborList {
public:
    using Index = util::Index<DIM>;
    NeighborList(std::array<dtype, DIM> gridSize, dtype interactionRadius, int nSubdivides = 2)
            : _gridSize(gridSize), nSubdivides(nSubdivides) {
        std::array<std::uint32_t, DIM> nCells;
        for (int i = 0; i < DIM; ++i) {
            _cellSize[i] = interactionRadius / nSubdivides;
            if (gridSize[i] <= 0) throw std::invalid_argument("grid sizes must be positive.");
            nCells[i] = gridSize[i] / _cellSize[i];
        }
        _index = Index(nCells);
        head.resize(_index.size());

        std::vector<std::size_t> adj;
        adj.reserve(1 + std::pow(nSubdivides, DIM));
        for(std::size_t i = 0; i < nCellsTotal(); ++i) {
            auto ijk = _index.inverse(i);

        }
    }

    ~NeighborList() = default;
    NeighborList(const NeighborList&) = delete;
    NeighborList &operator=(const NeighborList&) = delete;

    NeighborList(NeighborList&&) noexcept = default;
    NeighborList &operator=(NeighborList&&) noexcept = default;

    void setTypes(const std::unordered_set<std::size_t> &allowedTypes) {
        if constexpr(!allTypes) {
            types = allowedTypes;
        }
    }

    constexpr bool isAllowedType(const auto &type) const {
        if constexpr(allTypes) {
            return true;
        } else {
            return types.find(type) != end(types);
        }
    }

    template<typename ParticleCollection, typename Pool>
    void update(std::shared_ptr<ParticleCollection> collection, std::shared_ptr<Pool> pool) {
        list.resize(collection->size() + 1);
        std::fill(std::begin(list), std::end(list), 0);
        std::fill(std::begin(head), std::end(head), thread::copyable_atomic<std::size_t>());
        auto updateOp = [this](std::size_t particleId, const auto &pos, const auto &type, const auto&/*noe*/) {
            if(isAllowedType(type)) {
                auto boxId = positionToBoxIx(&pos.data[0]);

                // CAS
                auto &atomic = *head.at(boxId);
                auto currentHead = atomic.load();
                while (!atomic.compare_exchange_weak(currentHead, particleId)) {}
                list[particleId] = currentHead;
            }
        };

        collection->forEachParticle(updateOp, pool);
        pool->waitForTasks();
    }

    typename Index::GridDims gridPos(const dtype *pos) const {
        std::array<std::uint32_t, DIM> projections;
        for (auto i = 0u; i < DIM; ++i) {
            projections[i] = static_cast<std::uint32_t>(
                    std::max(static_cast<dtype>(0.), static_cast<dtype>(std::floor((pos[i] + .5 * _gridSize[i]) / _cellSize[i]))));
            projections[i] = util::clamp(projections[i], 0U, _index[i] - 1);
        }
        return projections;
    }

    std::uint32_t positionToBoxIx(const dtype *pos) const {
        auto boxId = _index.index(gridPos(pos));
        return boxId;
    }

    template<typename T>
    std::uint32_t cellNeighborsBegin(const T &cellIndex, auto dim) const {
        if constexpr(!periodic) {
            return cellIndex[dim] >= nSubdivides ? (cellIndex[dim] - nSubdivides) : 0;
        } else {
            return cellIndex[dim] >= nSubdivides ? cellIndex[dim] - nSubdivides : (_index[dim]  - cellIndex[dim]) % _index[dim];
        }
    }

    template<typename T>
    std::uint32_t cellNeighborsEnd(const T &cellIndex, auto dim) const {
        if constexpr(!periodic) {
            return cellIndex[dim] < _index[dim] - nSubdivides ? (cellIndex[dim] + nSubdivides + 1) : _index[dim];
        } else {
            return (cellIndex[dim] + nSubdivides + 1) % _index[dim];
        }
    }

    template<typename F, typename PoolPtr>
    void forEachCell(F &&f, PoolPtr pool) const {
        auto worker = [op = std::forward<F>(f)](auto cellIndex) {
            op(cellIndex);
        };
        for(decltype(nCellsTotal()) i = 0; i < nCellsTotal(); ++i) {
            pool->push(worker, _index.inverse(i));
        }
    }

    auto nCellsTotal() const {
        return _index.size();
    }

    template<typename F>
    void forEachNeighborInCell(F &&f, typename Index::value_type cellIndex) const {
        forEachNeighborInCell(std::forward<F>(f), _index.inverse(cellIndex));
    }

    template<bool all, typename F>
    void forEachNeighborInCell(F &&f, const typename Index::GridDims &cellIndex) const {
        auto particleId = (*head.at(_index.index(cellIndex))).load();
        while (particleId != 0) {

            for(std::size_t d = 0; d < DIM; ++d) {
                auto begin = cellNeighborsBegin(cellIndex, d);
                auto end = cellNeighborsEnd(cellIndex, d);

                auto cellPos = cellIndex;
                for (auto k = begin; k != end; k = periodic ? (k + 1) % _index[d] : k + 1) {
                    cellPos[d] = k;

                    auto neighborCellId = _index.index(cellPos);
                    auto neighborId = (*head.at(neighborCellId)).load();
                    while (neighborId != 0) {
                        if constexpr(all) {
                            if(neighborId != particleId) {
                                f(particleId, neighborId);
                            }
                        } else {
                            if (neighborId > particleId) {
                                f(particleId, neighborId);
                            }
                        }
                        neighborId = list.at(neighborId);
                    }
                }
            }

            particleId = list.at(particleId);
        }
    }



    template<typename ParticleCollection, typename F>
    void forEachNeighbor(std::size_t id, ParticleCollection &collection, F &&fun) const {
        const auto &pos = collection.position(id);
        auto gridPos = this->gridPos(&pos[0]);
        for (int i = 0u; i < DIM; ++i) {
            auto begin = cellNeighborsBegin(gridPos, i);
            auto end = cellNeighborsEnd(gridPos, i);

            typename Index::GridDims cellPos {gridPos};
            for (auto k = begin; k != end; k = periodic ? (k + 1) % _index[i] : k + 1) {
                cellPos[i] = k;

                if (cellPos != gridPos) {
                    auto cellId = _index.index(cellPos);
                    auto neighborId = (*head.at(cellId)).load();
                    while (neighborId != 0) {
                        fun(neighborId, collection.position(neighborId), collection.typeOf(neighborId),
                            collection.force(neighborId));
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
                    fun(neighborId, collection.position(neighborId), collection.typeOf(neighborId),
                        collection.force(neighborId));
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
    int nSubdivides;
    std::unordered_set<std::size_t> types {};
};

}
