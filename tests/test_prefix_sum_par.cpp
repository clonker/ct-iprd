//
// Created by mho on 4/14/22.
//

#include <catch2/catch.hpp>

#include <random>
#include <vector>

#include <ctiprd/config.h>
#include <ctiprd/cpu/PrefixSum.h>

TEST_CASE("Prefix sum", "[utils]") {
    auto pool = ctiprd::config::make_pool(8);

    std::vector<int> v(7);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(1, 60);
    std::generate(begin(v), end(v), [&]() {
        return distrib(gen);
    });

    std::fill(begin(v), end(v), 1);

    {
        auto vcopy = ctiprd::cpu::PrefixSum::compute(begin(v), end(v),  pool);

        auto vcopy2 = v;
        std::exclusive_scan(begin(vcopy2), end(vcopy2), begin(vcopy2), 0);
        REQUIRE(vcopy == vcopy2);
    }

}
