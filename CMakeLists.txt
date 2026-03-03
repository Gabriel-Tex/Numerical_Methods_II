cmake_minimum_required(VERSION 3.28)
project(NumericalMethods)

add_executable(methods main.cpp src/diff.cpp src/function.cpp)
target_include_directories(methods PRIVATE include)
