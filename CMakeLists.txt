cmake_minimum_required(VERSION 3.10)

project(yrx_hw1 C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

add_executable(yrx_hw1
    src/main.c
    src/algebra.c
)
include_directories(inc)