//
// Created by mho on 2/2/22.
//

#include <fmt/format.h>
#include <catch2/catch.hpp>
#include <ctiprd/thread/ctpl.h>
#include <ctiprd/thread/bsho.h>
#include <ctiprd/thread/ops.h>

TEMPLATE_TEST_CASE("Thread pool sanity", "[pool]", ctpl::thread_pool, bsho::thread_pool) {
    TestType pool {5};

    REQUIRE(pool.size() == 5);

    std::vector<int> vec(1000, 0);
    std::vector<int> vecRef(1000, 1);
    std::vector<std::future<void>> futures;
    futures.reserve(vec.size());
    for(std::size_t i = 0; i < vec.size(); ++i) {
        futures.push_back(pool.push([&vec, loopIndex = i] {
            vec[loopIndex] = 1;
            {
                using namespace std::chrono_literals;
                std::this_thread::sleep_for(1ms);
            }
        }));
    }
    REQUIRE(futures.size() == vec.size());
    std::for_each(begin(futures), end(futures), [](auto &future) {future.wait();});
    REQUIRE_THAT(vec, Catch::Matchers::UnorderedEquals(vecRef));
}
