cmake_minimum_required(VERSION 3.20)

project(RaMPC)

enable_testing()

add_compile_options(-Ofast -g)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(mex)
