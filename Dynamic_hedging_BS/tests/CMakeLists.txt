add_executable(hedging_tests test_bs_price.cpp)
target_link_libraries(hedging_tests PRIVATE hedging GTest::gtest_main project_warnings core)

include(GoogleTest)
gtest_discover_tests(hedging_tests)