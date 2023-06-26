cmake_minimum_required(VERSION 3.14)

project(tatami_r
    VERSION 1.0.0
    DESCRIPTION "tatami bindings for abstract R matrices"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(tatami_r INTERFACE)

add_subdirectory(extern)

target_link_libraries(tatami_r INTERFACE tatami)

# Switch between include directories depending on whether the downstream is
# using the build directly or is using the installed package.
target_include_directories(tatami_r
    INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
