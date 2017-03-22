cmake_minimum_required(VERSION 3.5)
project(ExtendedKF)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp
        FusionEKF.cpp
        kalman_filter.cpp
        main.cpp
        tools.cpp)

add_executable(ExtendedKF ${SOURCE_FILES})

file(COPY ${CMAKE_SOURCE_DIR}/sample-laser-radar-measurement-data-1.txt DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/sample-laser-radar-measurement-data-2.txt DESTINATION ${CMAKE_BINARY_DIR})