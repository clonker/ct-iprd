//
// Created by mho on 3/15/22.
//
#pragma once

#include <ctiprd/thread/bsho.h>

namespace ctiprd::config {

using ThreadPool = bsho::thread_pool;

template<typename Pool>
auto threadGranularity(const Pool &pool) {
    return 2 * pool.size();
}

template<typename Pool>
auto threadGranularity(std::shared_ptr<Pool> pool) {
    return 2 * pool->size();
}

}
