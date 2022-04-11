#include <catch2/catch.hpp>

#include <ctiprd/cpu/NeighborList.h>

TEST_CASE("Neighbor list sanity", "[neighbor_list]") {
    auto pool = std::make_shared<ctiprd::config::ThreadPool> (6);


    std::array<float, 2> gridSize {5, 5};
    ctiprd::cpu::nl::NeighborList<2, true, float> nl {gridSize, 1.};

}
