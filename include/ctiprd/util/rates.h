//
// Created by mho on 5/4/22.
//

#pragma once

#include <cmath>
#include <numbers>

namespace ctiprd::rates {

template<typename T>
T macroscopicRate(const T microRate, const T diff1, const T diff2, const T radius) {
    const auto kappa = std::sqrt(microRate / (diff1 + diff2));
    return static_cast<T>(4) * std::numbers::pi_v<T> * (diff1 + diff2) * radius
        * (static_cast<T>(1) - std::tanh(kappa * radius) / (kappa * radius));
}

}
