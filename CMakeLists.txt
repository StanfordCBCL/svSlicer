
cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-fopenmp")
PROJECT(svslicer)

find_package(VTK REQUIRED PATHS ~/vtk/build/lib64/cmake/vtk-9.2/)
include(${VTK_USE_FILE})

add_executable(svslicer src/main.cpp)
target_link_libraries(svslicer ${VTK_LIBRARIES})
