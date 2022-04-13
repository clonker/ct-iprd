//
// Created by mho on 4/11/22.
//

#pragma once

#include <cstddef>
#include <tuple>
#include <type_traits>

namespace ctiprd::util::hash {

namespace detail{
template<typename TupR, typename Tup = std::remove_reference_t<TupR>, auto N = std::tuple_size_v<Tup>>
constexpr auto ReverseTuple(TupR&& t) {
    return [&t]<auto... I>(std::index_sequence<I...>) {
        constexpr std::array indices{(N - 1 - I)...};
        return std::tuple<std::tuple_element_t<indices[I], Tup>...>{
                std::get<indices[I]>(std::forward<TupR>(t))...
        };
    }(std::make_index_sequence<N>{});
}
}

/**
 * Simplified version of boost hash combine.
 * @param seed the seed
 * @param v the value
 */
template<typename T>
void combine(std::size_t &seed, const T &v) {
    seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template<typename T, typename Tup = std::remove_reference_t<T>, auto N = std::tuple_size_v<Tup>>
struct ForwardBackwardTupleHasher {
    /**
     * Evaluates the hash of tuples independent of their reversedness.
     * @param tuple the tuple to compare
     * @return the hash of the reversed tuple if tuple[0] > tuple[-1], otherwise hash of the non-reversed version
     */
    std::size_t operator()(const T &tuple) const {
        std::size_t seed{0};
        [&seed, &tuple]<auto... I>(std::index_sequence<I...>) {
            if (std::get<0>(tuple) > std::get<N - 1>(tuple)) {
                constexpr std::array is{I...};
                (combine(seed, std::get<is[I]>(tuple)),...);
            } else {
                constexpr std::array is{(N - 1 - I)...};
                (combine(seed, std::get<is[I]>(tuple)),...);
            }
        }(std::make_index_sequence<N>{});
        return seed;
    }
};

template<typename Tuple>
struct TupleCompare {
    bool operator()(const Tuple &tup1, const Tuple &tup2) {
        return (std::get<0>(tup1) <= std::get<1>(tup1) ? tup1 : detail::ReverseTuple(tup1))
                < std::get<0>(tup2) <= std::get<1>(tup2) ? tup2 : detail::ReverseTuple(tup2);
    }
};

template<typename T>
class ForwardBackwardTupleEquality {
public:
    /**
     * Evaluates the equality of two tuples independently of their reversedness.
     * @param lhs the one tuple
     * @param rhs the other tuple
     * @return true if lhs == rhs or lhs == reverse(rhs)
     */
    constexpr bool operator()(const T &lhs, const T &rhs) const {
        return lhs == rhs || lhs == detail::ReverseTuple(rhs);
    }
};

}
