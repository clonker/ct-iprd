//
// Created by mho on 3/15/22.
//

#pragma once

namespace ctiprd {

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

};

}
