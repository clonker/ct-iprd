//
// Created by mho on 5/5/21.
//

#pragma once

#include <algorithm>
#include <atomic>
#include <thread>
#include <unordered_set>
#include <type_traits>

#include <ctiprd/config.h>
#include <ctiprd/util/ops.h>
#include <ctiprd/util/Index.h>
#include <ctiprd/thread/utils.h>

namespace ctiprd::cpu::nl {

template<typename Index, bool periodic>
struct CellAdjacency {
    static_assert(std::is_signed_v<typename Index::value_type>, "Needs index to be signed!");

    bool withinBounds(const Index &index, const typename Index::GridDims &ijk) const {
        const bool nonNegative = std::all_of(begin(ijk), end(ijk), [](const auto &element) { return element >= 0; });
        bool nonExceeding = true;
        for(std::size_t dim = 0; dim < Index::Dims; ++dim) {
            nonExceeding &= ijk[dim] < index[dim];
        }
        return nonNegative && nonExceeding;
    }

    void findEachAdjCell(const Index &index, const typename Index::GridDims &ijk, std::size_t dim, std::int32_t r,
                         std::vector<std::size_t> &adj) const {
        for (int ii = ijk[dim] - r; ii <= ijk[dim] + r; ++ii) {
            typename Index::GridDims cix {ijk};
            cix[dim] = ii;
            if constexpr(periodic) {
                cix[dim] = (cix[dim] % index[dim] + index[dim]) % index[dim];  // wrap around
            }

            if(withinBounds(index, cix)) {
                if (dim > 0) {
                    findEachAdjCell(index, cix, dim - 1, r, adj);
                } else {
                    adj.push_back(index.index(cix));
                }
            }
        }
    }
    CellAdjacency() = default;
    CellAdjacency(const Index &index, const std::int32_t radius) {
        typename Index::GridDims nNeighbors {};
        for (std::size_t i = 0; i < Index::Dims; ++i) {
            nNeighbors[i] = std::min(index[i], static_cast<decltype(index[i])>(2 * radius + 1));
        }
        // number of neighbors plus the cell itself
        const auto nAdjacentCells = 1 + std::accumulate(begin(nNeighbors), end(nNeighbors), 1, std::multiplies<>());
        // number of neighbors plus cell itself plus padding to save how many neighbors actually
        cellNeighbors = util::Index<2>{std::array<typename Index::value_type, 2>{index.size(), static_cast<typename Index::value_type>(1 + nAdjacentCells)}};
        cellNeighborsContent.resize(cellNeighbors.size());

        std::vector<std::size_t> adj;
        adj.reserve(1 + nAdjacentCells);
        for(std::size_t i = 0; i < index.size(); ++i) {
            adj.resize(0);
            const auto ijk = index.inverse(i);
            findEachAdjCell(index, ijk, Index::Dims - 1, radius, adj);
            std::sort(adj.begin(), adj.end());
            adj.erase(std::unique(std::begin(adj), std::end(adj)), std::end(adj));

            const auto begin = cellNeighbors(i, 0);
            cellNeighborsContent[begin] = adj.size();
            std::copy(adj.begin(), adj.end(), &cellNeighborsContent.at(begin + 1));
        }
    }

    CellAdjacency(const CellAdjacency&) = delete;
    CellAdjacency &operator=(const CellAdjacency&) = delete;

    CellAdjacency(CellAdjacency&&) noexcept = default;
    CellAdjacency &operator=(CellAdjacency&&) noexcept = default;

    ~CellAdjacency() = default;

    [[nodiscard]] auto nNeighbors(auto cellIndex) const {
        return cellNeighborsContent[cellNeighbors(cellIndex, 0)];
    }

    [[nodiscard]] auto cellsBegin(auto cellIndex) const {
        return begin(cellNeighborsContent) + cellNeighbors(cellIndex, 1);
    }

    [[nodiscard]] auto cellsEnd(auto cellIndex) const {
        return cellsBegin(cellIndex) + nNeighbors(cellIndex);
    }

    // index of size (n_cells x (1 + nAdjacentCells)), where the first element tells how many adj cells are stored
    util::Index<2> cellNeighbors {};
    // backing vector of _cellNeighbors index of size (n_cells x (1 + nAdjacentCells))
    std::vector<std::size_t> cellNeighborsContent {};
};

template<int DIM, bool periodic, typename dtype, bool allTypes=true>
class NeighborList {
public:
    using Index = util::Index<DIM, std::array<std::int32_t, DIM>>;
    NeighborList(std::array<dtype, DIM> gridSize, dtype interactionRadius, int nSubdivides = 2)
            : _gridSize(gridSize), nSubdivides(nSubdivides) {
        std::array<typename Index::value_type, DIM> nCells;
        for (int i = 0; i < DIM; ++i) {
            _cellSize[i] = interactionRadius / nSubdivides;
            if (gridSize[i] <= 0) {
                throw std::invalid_argument("grid sizes must be positive.");
            }
            nCells[i] = gridSize[i] / _cellSize[i];
        }
        _index = Index(nCells);
        head.resize(_index.size());
        _adjacency = CellAdjacency<Index, periodic>{_index, nSubdivides};
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
        auto updateOp = [this](std::size_t particleId, const auto &pos, const auto &type, const auto&/*velocity*/) {
            if(isAllowedType(type)) {
                const auto boxId = positionToBoxIx(&pos.data[0]);

                // CAS
                auto &atomic = *head.at(boxId);
                auto currentHead = atomic.load();
                while (!atomic.compare_exchange_weak(currentHead, particleId + 1)) {}
                list[particleId + 1] = currentHead;
            }
        };

        const auto futures = collection->forEachParticle(updateOp, pool);
        for(const auto &future : futures) {
            future.wait();
        }
    }

    template<typename Position>
    typename Index::GridDims gridPos(const Position &pos) const {
        typename Index::GridDims projections;
        for (auto i = 0U; i < DIM; ++i) {
            projections[i] = static_cast<typename Index::value_type>(
                    std::max(static_cast<dtype>(0.), static_cast<dtype>(std::floor((pos[i] + .5 * _gridSize[i]) / _cellSize[i]))));
            projections[i] = std::clamp(projections[i], static_cast<typename Index::value_type>(0),
                                        static_cast<typename Index::value_type>(_index[i] - 1));
        }
        return projections;
    }

    template<typename Position>
    std::uint32_t positionToBoxIx(const Position &pos) const {
        return _index.index(gridPos(pos));
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
    std::vector<std::future<void>> forEachCell(F &&func, PoolPtr pool) const {
        std::vector<std::future<void>> futures;
        const auto worker = [operation = std::forward<F>(func)](auto begin, auto end) {
            for(auto i = begin; i != end; ++i) {
                operation(i);
            }
        };
        const auto granularity = config::threadGranularity(pool);
        const auto grainSize = nCellsTotal() / granularity;
        std::size_t grainStep {};
        for (grainStep = 0; grainStep < granularity - 1; ++grainStep) {
            futures.emplace_back(pool->push(worker, grainStep * grainSize, (grainStep + 1) * grainSize));
        }
        if(grainStep * grainSize != nCellsTotal()) {
            futures.emplace_back(pool->push(worker, grainStep * grainSize, nCellsTotal()));
        }
        return futures;
    }

    auto nCellsTotal() const {
        return _index.size();
    }

    template<typename F>
    void forEachNeighborInCell(F &&func, const typename Index::GridDims &cellIndex) const {
        forEachNeighborInCell(std::forward<F>(func), _index(cellIndex));
    }

    template<bool all, typename F>
    void forEachNeighborInCell(F &&func, typename Index::value_type cellIndex) const {
        auto particleId = (*head.at(cellIndex)).load();
        while (particleId != 0) {
            const auto begin = _adjacency.cellsBegin(cellIndex);
            const auto end = _adjacency.cellsEnd(cellIndex);
            for (auto k = begin; k != end; ++k) {
                const auto neighborCellId = *k;
                auto neighborId = (*head.at(neighborCellId)).load();
                while (neighborId != 0) {
                    if constexpr(all) {
                        if(neighborId != particleId) {
                            func(particleId - 1, neighborId - 1);
                        }
                    } else {
                        if (neighborId > particleId) {
                            func(particleId - 1, neighborId - 1);
                        }
                    }
                    neighborId = list.at(neighborId);
                }
            }

            particleId = list.at(particleId);
        }
    }



    template<typename ParticleCollection, typename F>
    void forEachNeighbor(std::size_t particleId, ParticleCollection &collection, F &&fun) const {
        const auto &pos = collection.position(particleId);
        const auto gridPos = this->gridPos(&pos[0]);
        const auto gridId = _index.index(gridPos);
        const auto begin = _adjacency.cellsBegin(gridId);
        const auto end = _adjacency.cellsEnd(gridId);
        for (auto k = begin; k != end; ++k) {
            const auto neighborCellId = *k;

            auto neighborId = (*head.at(neighborCellId)).load();
            while (neighborId != 0) {
                if (neighborId - 1 != particleId) {
                    fun(neighborId - 1, collection.position(neighborId - 1), collection.typeOf(neighborId - 1),
                        collection.force(neighborId - 1));
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
    Index _index{};
    CellAdjacency<Index, periodic> _adjacency {};
    int nSubdivides;
    std::unordered_set<std::size_t> types {};
};

}
