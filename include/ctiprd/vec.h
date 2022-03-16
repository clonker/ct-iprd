//
// Created by mho on 3/15/22.
//

#pragma once

namespace ctiprd {

namespace detail {
template <typename T> concept arithmetic = std::is_arithmetic_v<T>;
}

template<typename dtype, int DIM>
struct Vec {
    using data_type = std::array<dtype, DIM>;
    using size_type = typename data_type::size_type;
    using reference = typename data_type::reference;
    using value_type = typename data_type::value_type;
    using const_reference = typename data_type::const_reference ;
    data_type data;

    reference operator[](size_type pos) {
        return data[pos];
    }

    const_reference operator[](size_type pos) const {
        return data[pos];
    }

    bool operator==(const Vec &rhs) const {
        return data == rhs.data;
    }

    Vec& operator+=(const Vec &other) {
        for(size_type i = 0; i < DIM; ++i) data[i] += other[i];
        return *this;
    }

    template<typename T>
    Vec& operator+=(T arg) requires detail::arithmetic<T> {
        for(size_type i = 0; i < DIM; ++i) data[i] += arg;
        return *this;
    }

    template<typename T>
    Vec& operator-=(T arg) requires detail::arithmetic<T> {
        for(size_type i = 0; i < DIM; ++i) data[i] -= arg;
        return *this;
    }

    template<typename T>
    Vec& operator*=(T arg) requires detail::arithmetic<T> {
        for(size_type i = 0; i < DIM; ++i) data[i] *= arg;
        return *this;
    }

    template<typename T>
    Vec& operator/=(T arg) requires detail::arithmetic<T> {
        for(size_type i = 0; i < DIM; ++i) data[i] /= arg;
        return *this;
    }

    template<typename T>
    friend Vec operator*(Vec lhs, T rhs) requires detail::arithmetic<T> {
        lhs *= rhs;
        return lhs;
    }

    template<typename T>
    friend Vec operator/(Vec lhs, T rhs) requires detail::arithmetic<T> {
        lhs /= rhs;
        return lhs;
    }

    friend Vec operator+(Vec lhs, const Vec &rhs) {
        lhs += rhs;
        return lhs;
    }

};

}
