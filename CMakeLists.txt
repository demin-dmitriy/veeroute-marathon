cmake_minimum_required(VERSION 3.12.1)

project(veeroute-marathon)

add_executable(solver src/main.cpp)
target_compile_options(solver PUBLIC -std=c++17 -Wall)
target_link_options(solver PUBLIC -static)
