project(lv_mod)

pybind11_add_module(lv_mod lotka_volterra.cpp)
target_link_libraries(lv_mod PUBLIC ct-iprd::ct-iprd)
