/**
 *
 *
 * @file pbc.h
 * @brief 
 * @author clonker
 * @date 4/8/22
 */
#pragma once

namespace ctiprd::util::pbc {

template<typename System, typename Position>
void wrapPBC(Position &pos) {
    if constexpr(System::periodic) {
        for (int d = 0; d < System::DIM; ++d) {
            while (pos[d] >= .5 * System::boxSize[d]) pos[d] -= System::boxSize[d];
            while (pos[d] < -.5 * System::boxSize[d]) pos[d] += System::boxSize[d];
        }
    }
}

template<typename System, typename Position>
auto shortestDifference(const Position &p1, const Position &p2) {
    auto diff = p2 - p1;
    wrapPBC<System>(diff);
    return diff;
}

template<typename System, typename Position>
auto dSquared(const Position &p1, const Position &p2) {
    auto diff = shortestDifference<System>(p1, p2);
    return diff * diff;
}
}
