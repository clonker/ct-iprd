//
// Created by mho on 4/11/22.
//
#include <catch2/catch.hpp>

#include <ctiprd/util/hash.h>

TEST_CASE("Forward Backward Map", "[tuples][hash]") {
    using Key = std::tuple<std::size_t, std::size_t>;
    using Hasher = ctiprd::util::hash::ForwardBackwardTupleHasher<Key>;
    using Eq = ctiprd::util::hash::ForwardBackwardTupleEquality<Key>;

    std::unordered_map<Key, std::string, Hasher, Eq> map;

    Key key {3, 4};
    map.insert({key, "test"});

    REQUIRE(map.find(std::make_tuple(3, 4)) != end(map));
    REQUIRE(map.find(std::make_tuple(4, 3)) != end(map));
    REQUIRE(map.find(std::make_tuple(4, 3))->second == "test");
}
