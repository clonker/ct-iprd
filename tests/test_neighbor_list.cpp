//
// Created by mho on 4/22/22.
//
#include <catch2/catch.hpp>
#include <ctiprd/cpu/NeighborList.h>

TEST_CASE("Cell adjacency full", "[nl]") {
    using Index = ctiprd::util::Index<2, std::array<std::int32_t, 2>>;
    Index ix {15, 15};
    ctiprd::cpu::nl::CellAdjacency<Index, true> adj {ix, 2};

    for(auto i = 0; i < ix.size(); ++i) {
        REQUIRE(adj.nNeighbors(i) == 25);
    }
}

TEST_CASE("Cell adjacency full 3d", "[nl]") {
    using Index = ctiprd::util::Index<3, std::array<std::int32_t, 3>>;
    Index ix {15, 15, 15};
    ctiprd::cpu::nl::CellAdjacency<Index, true> adj {ix, 2};

    for(auto i = 0; i < ix.size(); ++i) {
        REQUIRE(adj.nNeighbors(i) == 125);
    }
}

TEST_CASE("Cell adjacency thin", "[nl]") {
    using Index = ctiprd::util::Index<2, std::array<std::int32_t, 2>>;
    Index ix {15, 3};
    ctiprd::cpu::nl::CellAdjacency<Index, true> adj {ix, 2};
    for(auto i = 0; i < ix.size(); ++i) {
        REQUIRE(adj.nNeighbors(i) == 5*3);
    }
}

TEST_CASE("Cell adjacency nonperiodic", "[nl]") {
    using Index = ctiprd::util::Index<2, std::array<std::int32_t, 2>>;
    Index ix {15, 15};
    ctiprd::cpu::nl::CellAdjacency<Index, false> adj {ix, 2};

    REQUIRE(adj.nNeighbors(0) == 9);
    REQUIRE(adj.nNeighbors(ix(5, 5)) == 25);
}
