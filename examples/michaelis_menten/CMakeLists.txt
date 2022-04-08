project(mm_mod)

pybind11_add_module(mm_mod michaelis_menten.cpp)
target_link_libraries(mm_mod PUBLIC ct-iprd::ct-iprd)
