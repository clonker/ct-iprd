find_package(benchmark REQUIRED)

add_executable(bench_map bench_map.cpp)
target_link_libraries(bench_map ct-iprd::ct-iprd benchmark::benchmark)
