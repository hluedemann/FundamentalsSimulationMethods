cmake_minimum_required(VERSION 3.14)

project(exercise1
    LANGUAGES CXX)

# Dependencies

# Targets

add_executable(part1 "part1.cpp")

# Properties
target_compile_features(part1
    PRIVATE
        cxx_std_17)