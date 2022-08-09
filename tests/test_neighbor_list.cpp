//
// Created by mho on 4/22/22.
//
#include <catch2/catch.hpp>
#include <ctiprd/cpu/NeighborList.h>
#include <ctiprd/cpu/ParticleCollection.h>
#include <ctiprd/systems/double_well.h>

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

TEST_CASE("Test NeighborList for two particles", "[nl]") {
    using System = ctiprd::systems::DoubleWell<float>;
    using CollectionType = ctiprd::cpu::ParticleCollection<System, ctiprd::cpu::particles::positions, ctiprd::cpu::particles::forces>;
    auto pool = ctiprd::config::make_pool(8);

    ctiprd::cpu::nl::NeighborList<2, false, float, true> nl {{10., 10.}, 1., 1};
    CollectionType collection{};
    collection.addParticle({{0, 0}}, "A");
    collection.addParticle({{0, 0}}, "A");

    nl.update(&collection, pool);

    nl.forEachNeighborInCell<true>([](auto pId, auto nId) {
        spdlog::error("pid {}, nid {}", pId, nId);
    }, static_cast<std::uint32_t>(55));

    nl.forEachNeighbor(0, collection, [](auto nId, const auto&, const auto&, const auto&) {
        spdlog::error("0 got neighbor {}", nId);
    });

    nl.forEachNeighbor(1, collection, [](auto nId, const auto&, const auto&, const auto&) {
        spdlog::error("1 got neighbor {}", nId);
    });
}
