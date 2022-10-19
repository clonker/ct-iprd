project(lv3d_mod)

pybind11_add_module(${PROJECT_NAME} lotka_volterra_3d.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC ct-iprd::ct-iprd)
