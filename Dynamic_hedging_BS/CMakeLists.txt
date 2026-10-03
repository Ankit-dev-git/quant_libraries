add_library(hedging
    DataContainerCreation.cpp
    GridCreation.cpp
    delta_hedge.cpp
    generate_paths.cpp
    parallel_for.cpp
)
target_include_directories(hedging PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
find_package(Threads REQUIRED)
target_link_libraries(hedging PUBLIC Threads::Threads PRIVATE project_warnings core)

add_executable(dynamic_hedging Dynamic_hedging_BS.cpp)
target_link_libraries(dynamic_hedging PRIVATE hedging project_warnings core)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()