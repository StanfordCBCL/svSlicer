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
# are wired up. Never search for VTK without components: VTK 9 then loads all
# installed modules and their dependencies (e.g. MPI, which needs the C
# language). An older VTK fails the version check without loading any
# modules, and is then requested with its "vtk"-prefixed module names.
find_package(VTK 9 QUIET COMPONENTS
    CommonCore
    CommonDataModel
    FiltersCore
    IOXML)
if(NOT VTK_FOUND)
    find_package(VTK REQUIRED COMPONENTS
        vtkCommonCore
        vtkCommonDataModel
        vtkFiltersCore
        vtkIOXML)
    include(${VTK_USE_FILE})
endif()
message(STATUS "Using VTK ${VTK_VERSION}")

add_executable(svslicer src/main.cpp)
target_link_libraries(svslicer PRIVATE ${VTK_LIBRARIES} OpenMP::OpenMP_CXX)

option(SVSLICER_BUILD_TESTS "Build the svSlicer regression tests" ON)
if(SVSLICER_BUILD_TESTS)
    enable_testing()
    add_executable(svslicer_test_utils tests/svslicer_test_utils.cpp)
    target_link_libraries(svslicer_test_utils PRIVATE ${VTK_LIBRARIES})

    # Multithreaded slicing must give exactly the single-threaded result.
    # More threads than cores are used on purpose, to provoke interleaving.
    set(SVSLICER_TEST_THREADS 8 CACHE STRING
        "Number of OpenMP threads used by the regression test")
    set(SVSLICER_TEST_REPEATS 5 CACHE STRING
        "Number of multithreaded runs in the regression test")
    add_test(NAME omp_matches_single_thread
        COMMAND ${CMAKE_COMMAND}
            -DSVSLICER=$<TARGET_FILE:svslicer>
            -DTEST_UTILS=$<TARGET_FILE:svslicer_test_utils>
            -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/test_omp
            -DTHREADS=${SVSLICER_TEST_THREADS}
            -DREPEATS=${SVSLICER_TEST_REPEATS}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_omp_test.cmake)

    add_test(NAME velocity_components_and_thread_cap
        COMMAND ${CMAKE_COMMAND}
            -DSVSLICER=$<TARGET_FILE:svslicer>
            -DTEST_UTILS=$<TARGET_FILE:svslicer_test_utils>
            -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/test_edge_cases
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_edge_tests.cmake)
endif()

# VTK 9 uses autoinit modules; this registers the required factory objects.
if(NOT VTK_VERSION VERSION_LESS "8.90")
    vtk_module_autoinit(TARGETS svslicer MODULES ${VTK_LIBRARIES})
    if(SVSLICER_BUILD_TESTS)
        vtk_module_autoinit(TARGETS svslicer_test_utils
            MODULES ${VTK_LIBRARIES})
    endif()
endif()

install(TARGETS svslicer DESTINATION bin)
