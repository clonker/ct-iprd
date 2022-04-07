add_executable(2d_double_well 2d_double_well.cpp)
target_link_libraries(2d_double_well PUBLIC ct-iprd::ct-iprd ct-iprd::ct-iprd CONAN_PKG::spdlog CONAN_PKG::catch2)

pybind11_add_module(dw_mod 2d_double_well_module.cpp)
target_link_libraries(dw_mod PUBLIC ct-iprd::ct-iprd ct-iprd::ct-iprd CONAN_PKG::spdlog CONAN_PKG::catch2 CONAN_PKG::pybind11)
