cmake_minimum_required(VERSION 3.17)
project(ProjectProgram)

set(CMAKE_CXX_STANDARD 11)

add_executable(ProjectProgram main.cpp src/solver.cpp src/solver.h src/solver_inputs.h src/Constants.h src/particle.h src/prng_engine.h src/engine.h)