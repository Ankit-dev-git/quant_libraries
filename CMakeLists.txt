cmake_minimum_required(VERSION 3.25)
project(quant_libraries CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(project_warnings INTERFACE)
target_compile_options(project_warnings INTERFACE -Wall -Wextra -Wpedantic)

include(CTest)
if(BUILD_TESTING)
    include(FetchContent)
    FetchContent_Declare(googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
        SYSTEM DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    FetchContent_MakeAvailable(googletest)
endif()

add_subdirectory(core)
add_subdirectory(Dynamic_hedging_BS)
add_subdirectory(FFT_solver)