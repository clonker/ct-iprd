//
// Created by mho on 4/5/22.
//

#pragma once

#include <string_view>
#include <cstddef>
#include <array>
#include <cmath>
#include <algorithm>
#include <fmt/format.h>

namespace ctiprd::geometry {

template<std::size_t dim, typename T>
struct Box {

    using dtype = T;
    static constexpr std::string_view name = "Box";

    template<bool inclusion, typename Position>
    [[nodiscard]] bool contains(const Position &position) const {
        bool result = true;
        #pragma unroll
        for(std::size_t d = 0; d < dim; ++d) {
            result &= position[d] > v0[d] && position[d] < v1[d];
        }
        if constexpr(!inclusion) {
            result = !result;
        }
        return result;
    }

    template<bool inclusion, typename Position>
    [[nodiscard]] Position smallestDifference(const Position &position) const {
        if constexpr(inclusion) {
            Position difference {};
            #pragma unroll
            for(std::size_t d = 0; d < dim; ++d) {
                if (position[d] < v0[d]) {
                    difference[d] = position[d] - v0[d];
                } else if (position[d] > v1[d]) {
                    difference[d] = position[d] - v1[d];
                }
            }

            return difference;
        } else {
            if(contains<true>(position)) {
                Position difference {};
                // All components of position are within range of (v0[d], v1[d]]. Find minimum.
                std::array<dtype, 2*dim> diffs;
                #pragma unroll
                for(std::size_t d = 0; d < dim; ++d) {
                    diffs[2*d] = std::abs(position[d] - v0[d]);
                    diffs[2*d+1] = std::abs(position[d] - v1[d]);
                }
                auto it = std::min_element(begin(diffs), end(diffs));
                auto ix = std::distance(begin(diffs), it);
                auto vix = ix % 2;

                if (vix == 0) {
                    difference[ix / 2] = position[ix / 2] - v0[ix / 2];
                } else {
                    difference[ix / 2] = position[ix / 2] - v1[ix / 2];
                }

                return difference;
            } else {
                return {};
            }
        }
    }

    [[nodiscard]] std::string describe() const {
        return fmt::format("minimum vertex v0={} and maximum vertex v1={}", v0, v1);
    }

    // lower left vertex and upper right vertex
    std::array<dtype, dim> v0 {}, v1 {};
};

}
