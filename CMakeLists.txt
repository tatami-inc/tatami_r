cmake_minimum_required(VERSION 3.24)

project(tatami_r
    VERSION 2.0.4
    DESCRIPTION "tatami bindings for abstract R matrices"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(tatami_r INTERFACE)
add_library(tatami::tatami_r ALIAS tatami_r)

option(TATAMI_R_FETCH_EXTERN "Automatically fetch tatami_r's external dependencies." ON)
if(TATAMI_R_FETCH_EXTERN)
    add_subdirectory(extern)
else()
    find_package(tatami_tatami 4.0.0 CONFIG REQUIRED)
    find_package(tatami_tatami_chunked 2.0.0 CONFIG REQUIRED)
    find_package(tatami_manticore 1.0.0 CONFIG REQUIRED)
    find_package(ltla_sanisizer 0.1.0 CONFIG REQUIRED)
endif()

target_link_libraries(tatami_r INTERFACE tatami::tatami tatami::tatami_chunked tatami::manticore ltla::sanisizer)

# Switch between include directories depending on whether the downstream is
# using the build directly or is using the installed package.
include(GNUInstallDirs)
target_include_directories(tatami_r
    INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/tatami_tatami_r>"
)

# Installing for find_package.
include(CMakePackageConfigHelpers)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tatami_tatami_r)

install(TARGETS tatami_r
    EXPORT tatami_rTargets)

install(EXPORT tatami_rTargets
    FILE tatami_tatami_rTargets.cmake
    NAMESPACE tatami::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami_r)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatami_rConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami_r)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/tatami_tatami_rConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatami_rConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/tatami_tatami_rConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/tatami_tatami_r)
