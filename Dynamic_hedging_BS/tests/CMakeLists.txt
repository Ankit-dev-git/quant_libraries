include(FetchContent)
FetchContent_Declare(googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_MakeAvailable(googletest)

add_executable(hedging_tests test_bs_price.cpp)
target_link_libraries(hedging_tests PRIVATE hedging GTest::gtest_main project_warnings)

include(GoogleTest)
gtest_discover_tests(hedging_tests)