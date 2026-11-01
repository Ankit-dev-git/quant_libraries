add_library(heston_mc
    qe_scheme.cpp
)

target_include_directories(heston_mc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(heston_mc PRIVATE project_warnings PUBLIC core)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()