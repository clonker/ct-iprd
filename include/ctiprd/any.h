//
// Created by mho on 3/28/22.
//

#pragma once

#include <tuple>

namespace ctiprd {

template<typename UnaryFunction, typename ...Types>
constexpr bool TupleAnyOf(UnaryFunction&& p, std::tuple<Types...>&& t) {
    return std::apply(
            [&](auto&& ...xs) constexpr { return (p(std::forward<decltype(xs)>(xs)) || ...); },
            std::forward<std::tuple<Types...>>(t));
}

}
