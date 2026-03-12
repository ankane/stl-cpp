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

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(BUILD_TESTING "" ON)
    option(BUILD_EXAMPLE "" ON)

    if(BUILD_TESTING)
        add_executable(test test/main.cpp test/mstl_test.cpp test/stl_test.cpp)
        target_link_libraries(test PRIVATE stl::stl)
        if(NOT MSVC)
            target_compile_options(test PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Werror)
        endif()
    endif()

    if(BUILD_EXAMPLE)
        add_executable(example example/main.cpp)
        target_link_libraries(example PRIVATE stl::stl)
    endif()
endif()
