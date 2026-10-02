cmake_minimum_required(VERSION 3.17)
project(ProjectProgram CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_library(lhmc src/solver.cpp)
target_include_directories(lhmc PUBLIC src)

add_executable(ProjectProgram main.cpp)
target_link_libraries(ProjectProgram PRIVATE lhmc)

# Verification test: run with `ctest --test-dir build --output-on-failure`
enable_testing()
add_executable(test_langmuir tests/test_langmuir.cpp)
target_link_libraries(test_langmuir PRIVATE lhmc)
add_test(NAME langmuir_verification COMMAND test_langmuir)
