//
// Created by mho on 3/15/22.
//
#pragma once

#include <thread>
#include <random>
#include <memory>

#include <spdlog/spdlog.h>

#include <ctiprd/thread/bsho.h>
// #include <ctiprd/thread/ctpl.h>

namespace ctiprd::config {

using ThreadPool = bsho::thread_pool;
template<typename Pool = ThreadPool>
using PoolPtr = std::shared_ptr<Pool>;

using DefaultGenerator = std::mt19937;

template<typename Pool = ThreadPool>
static auto make_pool(int n) {
    return std::make_shared<Pool>(n);
}

template<typename Pool>
auto threadGranularity(PoolPtr<Pool> pool) {
    static auto res = 4 * pool->size();
    return res;
}

}
