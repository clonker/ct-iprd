//
// Created by mho on 4/19/22.
//
#include <benchmark/benchmark.h>
#include <ctiprd/util/Index.h>
static void IndexInverseFoldingExpression(benchmark::State& state) {
    ctiprd::util::Index<3> ix{{10, static_cast<unsigned int>(state.range()), 5}};
    for (auto _ : state) {
        benchmark::DoNotOptimize(ix.inverse(33));
    }
}

BENCHMARK(IndexInverseFoldingExpression)->DenseRange(5, 15, 1);

template<typename Index, typename GridDims = typename Index::GridDims>
GridDims inverse(std::size_t idx, const Index &index) {
    GridDims result {};
    if (index.size() > 0) {
        auto prefactor = index.size() / index[0];
        for (std::size_t d = 0; d < Index::Dims - 1; ++d) {
            auto x = std::floor(idx / prefactor);
            result[d] = x;
            idx -= x * prefactor;
            prefactor /= index[d + 1];
        }
        result[Index::Dims - 1] = idx;
    }
    return result;
}

static void IndexInverseLoop(benchmark::State& state) {
    ctiprd::util::Index<3> ix{{10, static_cast<unsigned int>(state.range()), 5}};
    for (auto _ : state) {
        benchmark::DoNotOptimize(inverse(33, ix));
    }
}

BENCHMARK(IndexInverseLoop)->DenseRange(5, 15, 1);

BENCHMARK_MAIN();
