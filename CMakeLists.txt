cmake_minimum_required(VERSION 3.10)  
project(MVM VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(MVM STATIC
    src/vmath.cpp
    include/vmath.hpp
)

target_include_directories(MVM
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)