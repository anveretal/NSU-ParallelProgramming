cmake_minimum_required(VERSION 3.22)
project(Lab1 VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(ct1 consistent_task_1.cpp)
add_executable(ct2 consistent_task_2.cpp)