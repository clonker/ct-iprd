#pragma once

#include <utility>
#include <array>
#include <numeric>
#include <cmath>

namespace ctiprd::util {
namespace detail {
template<typename... Ix>
struct ComputeIndex {
    template<typename Strides, typename Indices = std::make_index_sequence<sizeof...(Ix)>>
    static constexpr auto compute(const Strides &strides, Ix &&... ix) {
        std::tuple<Ix...> tup(std::forward<Ix>(ix)...);
        return compute(strides, tup, Indices{});
    }

    template<typename Arr, std::size_t... I>
    static constexpr auto compute(const Arr &strides, const std::tuple<Ix...> &tup, std::index_sequence<I...>) {
        return (0 + ... + (strides[I] * std::get<I>(tup)));
    }

    template<typename Arr, typename Arr2, std::size_t... I>
    static constexpr auto computeContainer(const Arr &strides, const Arr2 &tup, std::index_sequence<I...>) {
        return (0 + ... + (strides[I] * std::get<I>(tup)));
    }

};

template<typename GridDims, typename Strides, std::size_t N = std::tuple_size_v<Strides>>
constexpr GridDims indexInverse(const Strides &strides, std::size_t ix) {
    /*GridDims result {};
    [&]<auto... I>(std::index_sequence<I...>) {
        ([&]() {
            result[I] = std::floor(ix / std::get<I>(strides));
            ix -= result[I] * std::get<I>(strides);
        }(), ...);
    }(std::make_index_sequence<N>{});*/
    GridDims result {};
    for (std::size_t d = 0; d < N - 1; ++d) {
        auto x = std::floor(ix / strides[d]);
        result[d] = x;
        ix -= x * strides[d];
    }
    result[N - 1] = ix;
    return result;
}

}

template<std::size_t DIMS, typename T = std::array<std::uint32_t, DIMS>>
class Index {
    static_assert(DIMS > 0, "Dims has to be > 0");
public:
    using GridDims = T;
    static constexpr std::size_t Dims = DIMS;
    /**
     * The value type, inherited from GridDims::value_type
     */
    using value_type = typename GridDims::value_type;

    template<typename It>
    static auto make_index(It shapeBegin, It shapeEnd) {
        GridDims dims{};
        std::copy(shapeBegin, shapeEnd, begin(dims));
        auto n_elems = std::accumulate(begin(dims), end(dims), static_cast<value_type>(1),
                                       std::multiplies<value_type>());

        GridDims strides;
        strides[0] = n_elems / dims[0];
        for (std::size_t d = 0; d < Dims - 1; ++d) {
            strides[d + 1] = strides[d] / dims[d + 1];
        }

        return Index<Dims, GridDims>{dims, strides, n_elems};
    }

    template<typename Container = std::initializer_list<typename GridDims::value_type>>
    static auto make_index(const Container &container) {
        return make_index(begin(container), end(container));
    }

    /**
     * Constructs an empty index object of specified dimensionality. Not of much use, really.
     */
    Index() : _size(), _cum_size(), n_elems(0) {}

    Index(std::initializer_list<value_type> ilist) : Index(std::vector<value_type>{ilist}) {

    }

    template<typename Shape>
    explicit Index(const Shape &size)
            : _size(), n_elems(std::accumulate(begin(size), end(size), 1u, std::multiplies<value_type>())) {
        std::copy(begin(size), end(size), begin(_size));

        GridDims strides;
        strides[0] = n_elems / size[0];
        for (std::size_t d = 0; d < Dims - 1; ++d) {
            strides[d + 1] = strides[d] / size[d + 1];
        }
        _cum_size = std::move(strides);
    }

    /**
     * Constructs an index object with a number of size_t arguments that must coincide with the number of dimensions,
     * specifying the grid.
     * @tparam Args the argument types, must all be size_t
     * @param args the arguments
     */
    Index(GridDims size, GridDims strides, value_type nElems) : _size(std::move(size)), _cum_size(std::move(strides)),
                                                                n_elems(nElems) {}

    /**
     * the number of elements in this index, exactly the product of the grid dimensions
     * @return the number of elements
     */
    value_type size() const {
        return n_elems;
    }

    /**
     * Retrieve size of N-th axis
     * @tparam N the axis
     * @return size of N-th axis
     */
    template<int N>
    constexpr value_type get() const {
        return _size[N];
    }

    /**
     * retrieve size of N-th axis
     * @param N N
     * @return size of N-th axis
     */
    template<typename D>
    constexpr value_type operator[](D N) const {
        return _size[N];
    }

    /**
     * map Dims-dimensional index to 1D index
     * @tparam Ix the d-dimensional index template param type
     * @param ix the d-dimensional index
     * @return the 1D index
     */
    template<typename... Ix>
    constexpr value_type operator()(Ix &&... ix) const {
        static_assert(std::size_t(sizeof...(ix)) == Dims, "wrong input dim");
        return detail::ComputeIndex<Ix...>::compute(_cum_size, std::forward<Ix>(ix)...);
    }

    /**
     * map Dims-dimensional array to 1D index
     * @param indices the Dims-dimensional index
     * @return the 1D index
     */
    template<typename Arr, typename Indices = std::make_index_sequence<Dims>>
    constexpr value_type index(const Arr &indices) const {
        return detail::ComputeIndex<>::computeContainer(_cum_size, indices, Indices{});
    }

    /**
     * Inverse mapping 1D index to Dims-dimensional tuple
     * @param idx
     * @return
     */
    GridDims inverse(std::size_t idx) const {
        return detail::indexInverse<GridDims>(_cum_size, idx);
    }

private:
    GridDims _size;
    GridDims _cum_size;
    value_type n_elems;
};

}
