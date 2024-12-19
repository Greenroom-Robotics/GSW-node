cmake_minimum_required(VERSION 3.30)
project(gsw-node VERSION 0.1.0)

add_library(
        gsw-node SHARED
        src/binding.cpp
        src/BoundTeosBase.cpp
        src/BoundTeosIce.cpp
        src/BoundTeosSea.cpp
        ${CMAKE_JS_SRC}
)
target_include_directories(gsw-node PRIVATE include)
target_compile_features(gsw-node PUBLIC cxx_std_20)
target_link_libraries(gsw-node PRIVATE teos-cpp)


if (CMAKE_JS_VERSION)
    target_link_libraries(gsw-node PRIVATE ${CMAKE_JS_LIB})
    target_include_directories(gsw-node PUBLIC SYSTEM ${CMAKE_JS_INC})
    set_target_properties(gsw-node PROPERTIES PREFIX "" SUFFIX ".node")

    if (MSVC AND CMAKE_JS_NODELIB_DEF AND CMAKE_JS_NODELIB_TARGET)
        # Generate node.src
        execute_process(COMMAND ${CMAKE_AR} /def:${CMAKE_JS_NODELIB_DEF} /out:${CMAKE_JS_NODELIB_TARGET} ${CMAKE_STATIC_LINKER_FLAGS})
    endif ()
endif ()