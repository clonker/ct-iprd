add_executable(2d_double_well 2d_double_well.cpp)
target_link_libraries(2d_double_well PUBLIC ${CONAN_LIBS} ct-iprd::ct-iprd)
