add_library(core 
    normal.cpp
    BS_price.cpp)

target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(core PRIVATE project_warnings)