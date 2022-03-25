project(ct-iprd-tests)

set(TEST_SRC test_sanity.cpp test_thread_pool.cpp test_particle_collection.cpp test_system.cpp)
add_executable(tests catch_main.cpp ${TEST_SRC})
set_target_properties(tests PROPERTIES CXX_VISIBILITY_PRESET hidden)

target_link_libraries(tests PUBLIC ct-iprd::ct-iprd CONAN_PKG::spdlog CONAN_PKG::catch2)
if(OpenMP_FOUND)
    target_link_libraries(tests PUBLIC OpenMP::OpenMP_CXX)
endif()

add_custom_target(
        run_tests
        COMMAND tests --durations=yes -r console
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DEPENDS tests
)
