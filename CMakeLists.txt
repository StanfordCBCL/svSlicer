cmake_minimum_required(VERSION 3.16)
project(svSlicer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# On macOS the system AppleClang compiler does not ship the OpenMP runtime.
# Point CMake at the Homebrew "libomp" package so find_package(OpenMP) succeeds
# without having to hardcode a full LLVM toolchain.
if(APPLE)
    find_program(BREW_EXECUTABLE brew)
    if(BREW_EXECUTABLE)
        # Make Homebrew-installed packages (e.g. VTK) discoverable without the
        # caller having to pass -DCMAKE_PREFIX_PATH.
        execute_process(
            COMMAND ${BREW_EXECUTABLE} --prefix
            OUTPUT_VARIABLE HOMEBREW_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        if(HOMEBREW_PREFIX)
            list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}")
        endif()

        # VTK's CMake config lives under its own keg and is not symlinked into
        # the top-level Homebrew prefix, so add it explicitly.
        execute_process(
            COMMAND ${BREW_EXECUTABLE} --prefix vtk
            OUTPUT_VARIABLE VTK_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        if(VTK_PREFIX)
            list(APPEND CMAKE_PREFIX_PATH "${VTK_PREFIX}")
        endif()

        execute_process(
            COMMAND ${BREW_EXECUTABLE} --prefix libomp
            OUTPUT_VARIABLE LIBOMP_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        if(LIBOMP_PREFIX)
            set(OpenMP_ROOT "${LIBOMP_PREFIX}")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
            set(OpenMP_CXX_LIB_NAMES "omp")
            set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
        endif()
    endif()
endif()

find_package(OpenMP REQUIRED)

# VTK 9 requires the used modules to be requested explicitly, otherwise
# VTK_LIBRARIES is left empty and neither the headers nor the module autoinit
# are wired up.
find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    FiltersCore
    IOXML)

add_executable(svslicer src/main.cpp)
target_link_libraries(svslicer PRIVATE ${VTK_LIBRARIES} OpenMP::OpenMP_CXX)

# VTK 9 uses autoinit modules; this registers the required factory objects.
vtk_module_autoinit(TARGETS svslicer MODULES ${VTK_LIBRARIES})

install(TARGETS svslicer DESTINATION bin)
