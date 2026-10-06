add_library(recover_density
    density_recovery.cpp
    FFT.cpp
    characteristic_func.cpp
)

target_include_directories(recover_density PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(recover_density PRIVATE project_warnings)

add_executable(run_density_recovery run_density_recovery.cpp)
target_link_libraries(run_density_recovery PRIVATE recover_density project_warnings core)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()