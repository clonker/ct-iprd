project(lv_mod)

nanobind_add_module(lv_mod NOMINSIZE LTO lotka_volterra.cpp)
target_link_libraries(lv_mod PUBLIC ct-iprd::ct-iprd)
