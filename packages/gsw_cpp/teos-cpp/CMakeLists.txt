cmake_minimum_required(VERSION 3.30)
project(teos-cpp VERSION 0.1.0)

# Add the teos-cpp library
add_library(
        teos-cpp STATIC
        src/TeosBase.cpp
        src/TeosSea.cpp
        src/TeosIce.cpp
        src/SaarDataHandler.cpp
)
target_include_directories(teos-cpp PUBLIC include)
set_property(TARGET teos-cpp PROPERTY POSITION_INDEPENDENT_CODE ON)

add_subdirectory(tests/)
