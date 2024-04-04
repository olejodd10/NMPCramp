cmake_minimum_required(VERSION 3.20)

project(NMPCramp)

find_package(OpenMP)
message(STATUS "Adding \"${OpenMP_C_FLAGS}\" to CMAKE_C_FLAGS")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
add_compile_definitions(NUM_THREADS=8)

enable_testing()

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(mex)
