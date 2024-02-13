cmake_minimum_required(VERSION 3.26)
project(code_lab1)

set(CMAKE_CXX_STANDARD 17)

add_executable(code_lab1 main.cpp
        source.cpp
        method_QR.cpp
        method_gaussa.cpp
        header.h)
