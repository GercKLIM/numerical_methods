cmake_minimum_required(VERSION 3.26)
project(code_lab2)

set(CMAKE_CXX_STANDARD 17)

add_executable(code_lab2 main.cpp
        methods.cpp
        algebra.cpp
        algebra.h
        methods.h)
