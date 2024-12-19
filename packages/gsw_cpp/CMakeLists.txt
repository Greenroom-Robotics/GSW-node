cmake_minimum_required(VERSION 3.30)
project(gsw)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_definitions(-DNAPI_VERSION=9)

add_subdirectory(teos-cpp)
add_subdirectory(gsw-node)
