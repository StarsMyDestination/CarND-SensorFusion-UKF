cmake_minimum_required(VERSION 3.7)
project(UKF)

set(CMAKE_CXX_STANDARD 11)

set(sources
   src/ukf.cpp
   src/main.cpp
   src/tools.cpp)

add_executable(UKF ${sources})
