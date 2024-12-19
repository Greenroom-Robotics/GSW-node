Include(FetchContent)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.7.1
)

FetchContent_MakeAvailable(Catch2)

# Add a test excutable
add_executable(base src/base.cpp)
target_compile_features(base PUBLIC cxx_std_20)
target_link_libraries(base PRIVATE teos-cpp)
target_link_libraries(base PRIVATE Catch2::Catch2WithMain)

add_executable(sea src/sea.cpp)
target_compile_features(sea PUBLIC cxx_std_20)
target_link_libraries(sea PRIVATE teos-cpp)
target_link_libraries(sea PRIVATE Catch2::Catch2WithMain)

add_executable(ice src/ice.cpp)
target_compile_features(ice PUBLIC cxx_std_20)
target_link_libraries(ice PRIVATE teos-cpp)
target_link_libraries(ice PRIVATE Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(base REPORTER xml)
catch_discover_tests(sea REPORTER xml)
catch_discover_tests(ice REPORTER xml)