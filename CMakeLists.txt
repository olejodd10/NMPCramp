cmake_minimum_required(VERSION 3.20)

project(RaMPC)

enable_testing()

add_compile_options(-Ofast -g)

add_subdirectory(src)

find_package(Matlab)
if(Matlab_FOUND)
    add_subdirectory(mex)
endif()

add_subdirectory(test)
