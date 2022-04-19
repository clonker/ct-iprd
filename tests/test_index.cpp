//
// Created by mho on 4/19/22.
//
#include <catch2/catch.hpp>

#include <ctiprd/util/Index.h>

TEST_CASE("Inverse index", "[index]") {
    {
        ctiprd::util::Index<3> ix{{10, 3, 5}};
        REQUIRE(ix.inverse(ix(2, 1, 3))[0] == 2);
        REQUIRE(ix.inverse(ix(2, 1, 3))[1] == 1);
        REQUIRE(ix.inverse(ix(2, 1, 3))[2] == 3);
    }

    {
        ctiprd::util::Index<7> ix{{10, 3, 5, 7, 8, 2, 3}};
        REQUIRE(ix.inverse(ix(2, 1, 3, 2, 6, 1, 2))[0] == 2);
        REQUIRE(ix.inverse(ix(2, 1, 3, 2, 6, 1, 2))[1] == 1);
        REQUIRE(ix.inverse(ix(2, 1, 3, 2, 6, 1, 2))[2] == 3);
        REQUIRE(ix.inverse(ix(2, 1, 3, 2, 6, 1, 2))[3] == 2);
        REQUIRE(ix.inverse(ix(2, 1, 3, 2, 6, 1, 2))[4] == 6);
        REQUIRE(ix.inverse(ix(2, 1, 3, 2, 6, 1, 2))[5] == 1);
        REQUIRE(ix.inverse(ix(2, 1, 3, 2, 6, 1, 2))[6] == 2);
    }
}
