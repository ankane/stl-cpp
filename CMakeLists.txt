cmake_minimum_required(VERSION 3.18)

project(stl VERSION 0.2.0 LANGUAGES CXX)

include(GNUInstallDirs)

add_library(stl INTERFACE)
add_library(stl::stl ALIAS stl)

target_compile_features(stl INTERFACE cxx_std_20)

target_include_directories(
    stl
    INTERFACE
    "${PROJECT_SOURCE_DIR}/include/"
)

install(
    DIRECTORY "${PROJECT_SOURCE_DIR}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
