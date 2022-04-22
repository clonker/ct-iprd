//
// Created by mho on 2/1/22.
//

#pragma once

#include <type_traits>
#include <complex>
#include <cassert>
#include <numeric>

namespace ctiprd::util {

namespace detail {
template<typename T1, typename...>
struct variadic_first {
    /**
     * type of the first element of a variadic type tuple
     */
    using type = typename std::decay<T1>::type;
};

}

template<typename Iterator>
auto euclideanDistance(Iterator beginFirst, Iterator endFirst, Iterator beginSecond) {
    using dtype = typename std::iterator_traits<Iterator>::value_type;
    dtype init {0};
    for(; beginFirst != endFirst; ++beginFirst, ++beginSecond) {
        init += ((*beginFirst) - (*beginSecond)) * ((*beginFirst) - (*beginSecond));
    }
    return std::sqrt(init);
}

template<int dim, typename dtype>
dtype dot(const dtype *p1, const dtype *p2) {
    dtype result{0};
    for (auto i = 0u; i < dim; ++i) {
        result += p1[i] * p2[i];
    }
    return result;
}

template<int dim, typename dtype>
dtype length(const dtype *p) {
    dtype l{0};
    for (auto i = 0U; i < dim; ++i) {
        l += p[i] * p[i];
    }
    return std::sqrt(l);
}

template<typename T>
constexpr T factorial(T n) requires (std::is_integral_v<T>) {
    return n <= 1 ? 1 : factorial(n - 1) * n;
}

template<typename T, typename K>
constexpr T binomialCoefficient(T n, K k) requires (std::is_integral_v<T> && std::is_integral_v<K>) {
    return factorial(n) / (factorial(k) * factorial(n - k));
}

}
