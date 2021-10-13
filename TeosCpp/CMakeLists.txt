cmake_minimum_required(VERSION 3.0.0)
project(teos-cpp VERSION 0.1.0)

include(CTest)
enable_testing()

include_directories(
    include
)

# Add the teos-cpp library
add_library(teos-cpp STATIC
    src/TeosBase.cpp
    src/TeosSea.cpp
    src/TeosIce.cpp
    src/SaarDataHandler.cpp
)
set_property(TARGET teos-cpp PROPERTY POSITION_INDEPENDENT_CODE ON)

# Add a test excutable
add_executable(test-teos
    teos-cpp.cpp
)

# Link the library test code to lib
target_link_libraries(test-teos
    teos-cpp
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
