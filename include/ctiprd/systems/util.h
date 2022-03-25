//
// Created by mho on 3/25/22.
//

#pragma once

#include <string_view>
#include <algorithm>

namespace ctiprd::systems {

template<auto &types>
static constexpr std::size_t particleTypeId(std::string_view name) {
    std::size_t i = 0;
    for(auto it = begin(types); it != end(types); ++it, ++i) {
        if (it->name == name) {
            break;
        }
    }
    return i;
}

}
