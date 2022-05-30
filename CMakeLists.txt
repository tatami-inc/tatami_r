cmake_minimum_required(VERSION 3.14)

project(raticate
    VERSION 1.0.0
    DESCRIPTION "R bindings to the tatami library"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(raticate INTERFACE)

add_subdirectory(extern)

target_link_libraries(raticate INTERFACE tatami)
