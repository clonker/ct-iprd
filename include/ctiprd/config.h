//
// Created by mho on 3/15/22.
//
#pragma once

#include <ctiprd/thread/bsho.h>

namespace ctiprd::config {

using ThreadPool = bsho::thread_pool;
template<typename Pool = ThreadPool>
using PoolPtr = std::shared_ptr<Pool>;

template<typename Pool = ThreadPool>
static auto make_pool(int n) {
    return std::make_shared<Pool>(n);
}

template<typename Pool>
auto threadGranularity(PoolPtr<Pool> pool) {
    static auto res = 2 * pool->size();
    return res;
}

}
