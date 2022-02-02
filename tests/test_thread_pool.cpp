//
// Created by mho on 2/2/22.
//

#include <fmt/format.h>
#include <catch2/catch.hpp>
#include <ctiprd/thread/ctpl.h>

TEST_CASE("Thread pool sanity", "[ctpl]") {
    ctpl::ThreadPool pool {5};

    REQUIRE(pool.size() == 5);

    std::vector<int> vec(10000, 0);
    std::vector<int> vecRef(10000, 1);
    std::vector<std::future<void>> futures;
    futures.reserve(vec.size());
    for(std::size_t i = 0; i < vec.size(); ++i) {
        futures.push_back(pool.push([&vec, loopIndex = i](int threadIdx) {
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

