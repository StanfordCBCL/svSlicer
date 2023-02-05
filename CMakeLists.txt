
cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)
PROJECT(svSlicer)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-fopenmp, COMPILER_SUPPORTS_OPENMP)
if(COMPILER_SUPPORTS_OPENMP)
    add_compile_options(-fopenmp)
endif()

if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/qt5/")
endif()

find_package(VTK REQUIRED)
include_directories(${VTK_INCLUDE_DIRS})

add_executable(svslicer src/main.cpp)
target_link_libraries(svslicer ${VTK_LIBRARIES})
