
cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)
PROJECT(svSlicer)

if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/qt5/")
    set(CMAKE_CXX_COMPILER "/usr/local/opt/llvm/bin/clang")
    set(CMAKE_CXX_FLAGS "-I/usr/local/opt/llvm/include -fopenmp")
    set(CMAKE_EXE_LINKER_FLAGS "-L/usr/local/opt/llvm/lib")
else()
    set(CMAKE_CXX_FLAGS "-fopenmp")
endif()

find_package(VTK COMPONENTS IOXML FiltersCore REQUIRED)
include_directories(${VTK_INCLUDE_DIRS})

add_executable(svslicer src/main.cpp)
target_link_libraries(svslicer ${VTK_LIBRARIES})
install(TARGETS svslicer DESTINATION bin)
