# OPTIL.io uses CMake version 3.5.1, as of April 11, 2021.
# see https://www.optil.io/optilion/environment
cmake_minimum_required(VERSION 3.5)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# project
project(Exact)
add_executable(Exact src/exact.cpp)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)

# variables
add_compile_definitions(RANDOM_SEED=12345)
