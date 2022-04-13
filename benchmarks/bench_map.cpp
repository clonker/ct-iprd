#include <tuple>
#include <type_traits>
#include <map>
#include <unordered_map>
#include <random>

#include <benchmark/benchmark.h>

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

template<typename T>
void combine(std::size_t &seed, const T &v) {
    seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template<typename T, typename Tup = std::remove_reference_t<T>, auto N = std::tuple_size_v<Tup>>
class ForwardBackwardTupleHasher {
public:
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

template<typename T>
class ForwardBackwardTupleEquality {
public:
    constexpr bool operator()(const T &lhs, const T &rhs) const {
        return lhs == rhs || lhs == detail::ReverseTuple(rhs);
    }
};

template<typename Tuple>
struct TupleCompare {
    bool operator()(const Tuple &tup1, const Tuple &tup2) const {
        return (std::get<0>(tup1) <= std::get<1>(tup1) ? tup1 : detail::ReverseTuple(tup1))
               < (std::get<0>(tup2) <= std::get<1>(tup2) ? tup2 : detail::ReverseTuple(tup2));
    }
};

static void SortedMap(benchmark::State& state) {
    // Code inside this loop is measured repeatedly
    using Key = std::tuple<const int, const int>;
    std::map<Key, std::string, TupleCompare<Key>> map;
    for (int i = 0; i < state.range(); ++i) {
        for (int j = i+1; j < state.range(); ++j) {
            map[std::make_tuple(i, j)] = "huhu" + std::to_string(i) + std::to_string(j);
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, state.range()-1);

    for (auto _ : state) {
        state.PauseTiming();
        auto i = distrib(gen);
        auto j = distrib(gen);
        state.ResumeTiming();
        benchmark::DoNotOptimize(map[std::make_tuple(i, j)]);
    }
}
// Register the function as a benchmark
BENCHMARK(SortedMap)->DenseRange(1, 50, 1);

static void HashMap(benchmark::State& state) {
    // Code before the loop is not measured
    using Key = std::tuple<int, int>;
    std::unordered_map<Key, std::string, ForwardBackwardTupleHasher<Key>, ForwardBackwardTupleEquality<Key>> map;
    for (int i = 0; i < state.range(); ++i) {
        for (int j = i+1; j < state.range(); ++j) {
            map[std::make_tuple(i, j)] = "huhu" + std::to_string(i) + std::to_string(j);
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, state.range()-1);

    for (auto _ : state) {
        state.PauseTiming();
        auto i = distrib(gen);
        auto j = distrib(gen);
        state.ResumeTiming();
        benchmark::DoNotOptimize(map[std::make_tuple(i, j)]);
    }
}
BENCHMARK(HashMap)->DenseRange(1, 50, 1);


BENCHMARK_MAIN();
