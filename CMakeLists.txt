
cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fopenmp")
PROJECT(slice-results)

find_package(VTK REQUIRED PATHS ~/vtk/build/lib64/cmake/vtk-9.2/)
include(${VTK_USE_FILE})
set(VTK_SMP_IMPLEMENTATION_TYPE "OpenMP")
set(VTK_SMP_ENABLE_OPENMP ON)

add_executable(slicer main.cpp)
target_link_libraries(slicer ${VTK_LIBRARIES})
