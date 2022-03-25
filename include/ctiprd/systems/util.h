//
// Created by mho on 3/25/22.
//

#pragma once

#include <string_view>
#include <algorithm>

namespace ctiprd::systems {

template<const char* name, typename ParticleTypes>
struct ParticleTypeId {

};

template<auto types>
static constexpr std::size_t particleTypeId(std::string_view name) {
    std::size_t i = 0;
    std::for_each(begin(types), end(types), [name] (auto t){
        if constexpr(t.name == name) {

        }
    });
    return 0;
}

}
