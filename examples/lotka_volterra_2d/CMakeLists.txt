project(lv2d_mod)

pybind11_add_module(${PROJECT_NAME} lotka_volterra_2d.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC ct-iprd::ct-iprd)
