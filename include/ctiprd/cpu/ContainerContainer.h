/**
 *
 *
 * @file ContainerContainer.h
 * @brief 
 * @author clonker
 * @date 4/8/22
 */
#pragma once

namespace ctiprd::cpu {

template<typename... Containers>
struct ContainerContainerIterator {
    static_assert(sizeof...(Containers) > 0, "must not be empty");
    using InnerIterators = std::tuple<typename Containers::iterator...>;

    using iterator_category = std::random_access_iterator_tag;
    using difference_type = typename std::tuple_element_t<0, InnerIterators>::difference_type;
    using value_type = std::tuple<typename Containers::value_type...>;
    using pointer = value_type*;
    using reference = std::tuple<typename Containers::iterator::reference...>;

    InnerIterators innerIterators;

    bool operator==(const ContainerContainerIterator& r) const { return innerIterators == r.innerIterators; }
    bool operator!=(const ContainerContainerIterator& r) const { return innerIterators != r.innerIterators; }

    ContainerContainerIterator operator+(difference_type i) const {
        ContainerContainerIterator result {innerIterators};
        std::apply([&result, &i](auto &&... it) { ((it += i), ...); }, result.innerIterators);
        return result;
    }

    ContainerContainerIterator operator-(difference_type i) const {
        ContainerContainerIterator result {innerIterators};
        std::apply([&result, &i](auto &&... it) { ((it -= i), ...); }, result.innerIterators);
        return result;
    }

    difference_type operator-(const ContainerContainerIterator& r) const {
        return std::get<0>(innerIterators) - std::get<0>(r.innerIterators);
    }

    ContainerContainerIterator& operator++() {
        std::apply([](auto &&... it) { ((it++), ...); }, innerIterators);
        return *this;
    }

    ContainerContainerIterator& operator--() {
        std::apply([](auto &&... it) { ((it--), ...); }, innerIterators);
        return *this;
    }

    bool operator<(const ContainerContainerIterator& r) const {
        return std::get<0>(innerIterators) < std::get<0>(r.innerIterators);
    }

    reference operator*() {
        return [this]<auto... I>(std::index_sequence<I...>) {
            return reference{
                (*std::get<I>(innerIterators))...
            };
        } (std::make_index_sequence<sizeof...(Containers)>{});
    }
};

}
