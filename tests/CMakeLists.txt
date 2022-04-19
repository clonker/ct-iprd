project(ct-iprd-tests)

set(TEST_SRC
        test_sanity.cpp
        test_thread_pool.cpp
        test_particle_collection.cpp
        test_system.cpp
        test_reactions.cpp
        test_integration.cpp
        test_reactions_map.cpp
        test_cuda.cu
        test_michaelis_menten.cpp
        test_forward_backward_map.cpp
        test_prefix_sum_par.cpp
        test_index.cpp)
add_executable(tests catch_main.cpp ${TEST_SRC})
set_target_properties(tests PROPERTIES CXX_VISIBILITY_PRESET hidden)

target_link_libraries(tests PUBLIC ct-iprd::ct-iprd Catch2::Catch2)

add_custom_target(
        run_tests
        COMMAND tests --durations=yes -r console
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DEPENDS tests
)
