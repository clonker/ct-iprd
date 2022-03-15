//
// Created by mho on 3/15/22.
//

#pragma once

#include <type_traits>
#include <cstdint>
#include <atomic>
#include <vector>
#include <future>

namespace ctiprd::thread::ops {

template<typename T1, typename T2, typename F, typename Pool>
void parallelize_loop(const T1 &begin, const T2 &end, const F &loop, Pool &pool, std::uint32_t num_blocks = 0) {
    typedef std::common_type_t<T1, T2> T;
    T the_first_index = (T) begin;
    T last_index = (T) end;
    if (the_first_index == last_index)
        return;
    if (last_index < the_first_index) {
        T temp = last_index;
        last_index = the_first_index;
        the_first_index = temp;
    }
    last_index--;
    if (num_blocks == 0)
        num_blocks = pool.size();
    auto total_size = (std::size_t) (last_index - the_first_index + 1);
    auto block_size = (std::size_t) (total_size / num_blocks);
    if (block_size == 0) {
        block_size = 1;
        num_blocks = (std::uint32_t) total_size > 1 ? (std::uint32_t) total_size : 1;
    }
    std::vector<std::future<void>> futures;
    futures.reserve(num_blocks);
    for (std::uint32_t t = 0; t < num_blocks; t++) {
        T start = ((T) (t * block_size) + the_first_index);
        T end = (t == num_blocks - 1) ? last_index + 1 : ((T) ((t + 1) * block_size) + the_first_index);
        futures.push_back(pool.push([start, end, &loop] {
            loop(start, end);
        }));
    }
    std::for_each(begin(futures), end(futures), [](auto &future) {future.wait();});
}

}
