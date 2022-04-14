//
// Created by mho on 4/14/22.
//

#pragma once

#include <ctiprd/config.h>
#include <spdlog/spdlog.h>

namespace ctiprd::cpu {

struct PrefixSum {

    template<std::random_access_iterator InputIt, typename Pool>
    static auto compute(InputIt first, InputIt last, ctiprd::config::PoolPtr<Pool> pool) {
        std::vector out (first, last);
        /*auto granularity = ctiprd::config::threadGranularity(pool);
        auto n = std::distance(first, last);

        int offset = 1;
        for(int d = n >> 1; d > 0; d >>= 1) {
            for(int thid = 0; thid < d; ++thid) {
                int ai = offset * (2 * thid + 1) - 1;
                int bi = offset * (2 * thid + 2) - 1;

                out[bi] += out[ai];
            }
            offset *= 2;
        }

        out[n - 1] = 0;

        for(int d = 1; d < n; d *= 2) {
            offset >>= 1;

            for(int thid = 0; thid < d; ++thid) {
                int ai = offset*(2*thid+1)-1;
                int bi = offset*(2*thid+2)-1;

                float t = out[ai];
                out[ai] = out[bi];
                out[bi] += t;
            }
        }*/
        std::exclusive_scan(first, last, begin(out), 0);
        return out;
    }

};

}
