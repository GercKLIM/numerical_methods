cmake_minimum_required(VERSION 3.27)
project(Code_lab_1)

set(CMAKE_CXX_STANDARD 17)

add_executable(Code_lab_1 main.cpp
        algebra.h
        algebra.cpp
        methods.h
        methods.cpp
        tests.cpp)
