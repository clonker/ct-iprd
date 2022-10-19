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

/**
 * Struct computing and storing the neighboring cells in a CLL.
 *
 * @tparam Index Type of index that is used in the CLL.
 * @tparam periodic Whether periodic boundary conditions apply.
 */
template<typename Index, bool periodic>
struct CellAdjacency {
    static_assert(std::is_signed_v<typename Index::value_type>, "Needs index to be signed!");

    /**
     * Check whether a multi-index ijk is within bounds of the index itself.
     *
     * @param index The index.
     * @param ijk A multi index for which it is checked, whether it is still within bounds of index.
     * @return true if none of the index dimensions are exceeded and none of the entries are negative, else false
     */
    bool withinBounds(const Index &index, const typename Index::GridDims &ijk) const {
        const bool nonNegative = std::all_of(begin(ijk), end(ijk), [](const auto &element) { return element >= 0; });
        bool nonExceeding = true;
        for (std::size_t dim = 0; dim < Index::Dims; ++dim) {
            nonExceeding &= ijk[dim] < index[dim];
        }
        return nonNegative && nonExceeding;
    }

    /**
     * Finds adjacent cell indices given a current cell (specified by ijk) in a certain dimensional direction and
     * stores those inside an adjacency vector.
     *
     * @param index The reference index, specifying the grid.
     * @param ijk The current cell.
     * @param dim Dimension axis to explore.
     * @param r Number of subdivisions.
     * @param adj Vector with adjacent cells.
     */
    void findEachAdjCell(const Index &index, const typename Index::GridDims &ijk, std::size_t dim, std::int32_t r,
                         std::vector<std::size_t> &adj) const {
        for (int ii = ijk[dim] - r; ii <= ijk[dim] + r; ++ii) {
            typename Index::GridDims cix{ijk};
            cix[dim] = ii;
            if constexpr (periodic) {
                cix[dim] = (cix[dim] % index[dim] + index[dim]) % index[dim];  // wrap around
            }

            if (withinBounds(index, cix)) {
                if (dim > 0) {
                    findEachAdjCell(index, cix, dim - 1, r, adj);
                } else {
                    adj.push_back(index.index(cix));
                }
            }
        }
    }

    CellAdjacency() = default;

    /**
     * Creates and populates new cell adjacency struct.
     *
     * @param index The grid index object.
     * @param radius Radius in which to check for neighbors (in terms of discrete steps around a reference cell).
     */
    CellAdjacency(const Index &index, const std::int32_t radius) {
        typename Index::GridDims nNeighbors{};
        for (std::size_t i = 0; i < Index::Dims; ++i) {
            nNeighbors[i] = std::min(index[i], static_cast<decltype(index[i])>(2 * radius + 1));
        }
        // number of neighbors plus the cell itself
        const auto nAdjacentCells = 1 + std::accumulate(begin(nNeighbors), end(nNeighbors), 1, std::multiplies<>());
        // number of neighbors plus cell itself plus padding to save how many neighbors actually
        cellNeighbors = util::Index<2>{std::array<typename Index::value_type, 2>{
            index.size(), static_cast<typename Index::value_type>(1 + nAdjacentCells)}
        };
        cellNeighborsContent.resize(cellNeighbors.size());

        std::vector<std::size_t> adj;
        adj.reserve(1 + nAdjacentCells);
        for (std::size_t i = 0; i < index.size(); ++i) {
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

    CellAdjacency(const CellAdjacency &) = delete;

    CellAdjacency &operator=(const CellAdjacency &) = delete;

    CellAdjacency(CellAdjacency &&) noexcept = default;

    CellAdjacency &operator=(CellAdjacency &&) noexcept = default;

    ~CellAdjacency() = default;

    /**
     * The number of neighboring cells for specific (flat) cell index.
     *
     * @param cellIndex The cell index.
     * @return Number of neighbors.
     */
    [[nodiscard]] auto nNeighbors(auto cellIndex) const {
        return cellNeighborsContent[cellNeighbors(cellIndex, 0)];
    }

    /**
     * Random-access iterator denoting the begin of cell neighbors of cellIndex's cell.
     *
     * @param cellIndex The reference cell
     * @return iterator with begin
     */
    [[nodiscard]] auto cellsBegin(auto cellIndex) const {
        return begin(cellNeighborsContent) + cellNeighbors(cellIndex, 1);
    }

    /**
     * See cellsBegin.
     */
    [[nodiscard]] auto cellsEnd(auto cellIndex) const {
        return cellsBegin(cellIndex) + nNeighbors(cellIndex);
    }

    // index of size (n_cells x (1 + nAdjacentCells)), where the first element tells how many adj cells are stored
    util::Index<2> cellNeighbors{};
    // backing vector of _cellNeighbors index of size (n_cells x (1 + nAdjacentCells))
    std::vector<std::size_t> cellNeighborsContent{};
};

/**
 * A cell linked-list.
 *
 * @tparam DIM the dimension
 * @tparam periodic whether periodic boundary conditions apply
 * @tparam dtype the data type of positions (typically float/double)
 * @tparam allTypes whether this NL is for all types
 */
template<int DIM, bool periodic, typename dtype, bool allTypes=true>
class NeighborList {
public:
    /**
     * Index type which is used to access cells as spatial (i, j, k,...) indices or flat (ravel/unravel).
     */
    using Index = util::Index<DIM, std::array<std::int32_t, DIM>>;
    /**
     * Creates a new CLL based on a grid size (which assumed to result in an origin-centered space), an interaction
     * radius which is used to determine the amount of subdivision and a number of subdivides to further fine-grain
     * the sphere approximation.
     *
     * @param gridSize array of dtype describing the extent of space
     * @param interactionRadius maximum radius under which particles can interact
     * @param nSubdivides amount of fine-graining
     */
    NeighborList(std::array<dtype, DIM> gridSize, dtype interactionRadius, int nSubdivides = 2)
        : _gridSize(gridSize) {
            // determine the number of cells per axis
            std::array<typename Index::value_type, DIM> nCells;
            for (int i = 0; i < DIM; ++i) {
                _cellSize[i] = interactionRadius / nSubdivides;
                if (gridSize[i] <= 0) {
                    throw std::invalid_argument("grid sizes must be positive.");
                }
                nCells[i] = gridSize[i] / _cellSize[i];
            }
            // create index for raveling and unravling operations
            _index = Index(nCells);
            // initialize head to reflect the total number of cells
            head.resize(_index.size());
            // compute adjacencies among cells and store in arrays
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

    /**
     * Clears the currently stored particle indices and (re)computes the cell linked-list structure.
     *
     * @tparam Iterator Type of random-access iterator.
     * @param begin begin of the data structure containing vectors
     * @param end end of the data structure containing vectors
     * @param nJobs number of processes to use, threads are automatically joined
     */
    template<typename ParticleCollection, typename Pool>
    void update(ParticleCollection* collection, std::shared_ptr<Pool> pool) {
        // add artificial empty particle so that all entries can be unsigned
        list.resize(collection->size() + 1);
        // reset list
        std::fill(std::begin(list), std::end(list), 0);
        // reset head
        std::fill(std::begin(head), std::end(head), thread::copyable_atomic<std::size_t>());

        // operation which updates the cell linked-list by one individual particle based on its position
        const auto updateOp = [this](std::size_t particleId, const auto &pos, const auto &type, const auto&/*velocity*/) {
            if(isAllowedType(type)) {
                const auto boxId = positionToCellIndex(&pos.data[0]);

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

    /**
     * Determine a multidimensional cell index based on a position. This function is where the "space centered around
     * origin" assumption enters.
     *
     * @tparam Position Type of position
     * @param pos the position
     * @return multidimensional index
     */
    template<typename Position>
    typename Index::GridDims gridPos(const Position &pos) const {
        typename Index::GridDims projections;
        for (auto i = 0U; i < DIM; ++i) {
            projections[i] = static_cast<typename Index::value_type>(
                    std::max(static_cast<dtype>(0.),
                             static_cast<dtype>(std::floor((pos[i] + .5 * _gridSize[i]) / _cellSize[i]))));
            projections[i] = std::clamp(projections[i], static_cast<typename Index::value_type>(0),
                                        static_cast<typename Index::value_type>(_index[i] - 1));
        }
        return projections;
    }

    /**
     * Returns flat cell index based on position. See gridPos.
     *
     * @tparam Position Type of position.
     * @param pos The position
     * @return flat index pointing to a cell
     */
    template<typename Position>
    std::uint32_t positionToCellIndex(const Position &pos) const {
        return _index.index(gridPos(pos));
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

    /**
     * The number of cells in this cell linked-list
     */
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
    std::unordered_set<std::size_t> types {};
};

}
