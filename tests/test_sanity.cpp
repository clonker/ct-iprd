#include <catch2/catch.hpp>

#include <ctiprd/NeighborList.h>

TEST_CASE("Neighbor list sanity", "[neighbor_list]") {
    ctiprd::nl::NeighborList<2, true, float> nl;
}
