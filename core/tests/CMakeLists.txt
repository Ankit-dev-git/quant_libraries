add_executable(core_tests test_normal.cpp)
target_link_libraries(core_tests PRIVATE core GTest::gtest_main project_warnings)

include(GoogleTest)
gtest_discover_tests(core_tests)