cmake_minimum_required(VERSION 3.28)
project(NumericalMethods)

add_executable(methods main.cpp src/numericalDiff/finiteDiff.cpp 
src/numericalDiff/finiteDiffTaylor.cpp src/numericalDiff/finiteDiffNewton.cpp 
src/function.cpp)

target_include_directories(methods PRIVATE include)
