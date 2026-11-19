# Regression test for multithreaded slicing.
#
# Runs svSlicer on a synthetic tube once single-threaded and repeatedly with
# OMP_NUM_THREADS threads. Every output must match the analytic solution and
# agree with the single-threaded output to 1e-12 relative (summation-order
# rounding only). Concurrent use of a shared VTK dataset by the threads
# corrupts the results by 1e-3 and more, or crashes.
#
# Required variables: SVSLICER, TEST_UTILS, WORK_DIR, THREADS, REPEATS

foreach(var SVSLICER TEST_UTILS WORK_DIR THREADS REPEATS)
    if(NOT DEFINED ${var})
        message(FATAL_ERROR "${var} is not set")
    endif()
endforeach()

file(MAKE_DIRECTORY "${WORK_DIR}")
set(RESULT "${WORK_DIR}/tube.vtu")
set(CENTERLINE "${WORK_DIR}/tube_centerline.vtp")

function(run_checked)
    execute_process(COMMAND ${ARGN}
        WORKING_DIRECTORY "${WORK_DIR}"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE output)
    if(NOT result EQUAL 0)
        message("${output}")
        message(FATAL_ERROR "Command failed (${result}): ${ARGN}")
    endif()
    set(COMMAND_OUTPUT "${output}" PARENT_SCOPE)
endfunction()

function(run_svslicer threads output_file)
    set(ENV{OMP_NUM_THREADS} ${threads})
    run_checked("${SVSLICER}" "${RESULT}" "${CENTERLINE}" "${output_file}")
    set(SVSLICER_OUTPUT "${COMMAND_OUTPUT}" PARENT_SCOPE)
    run_checked("${TEST_UTILS}" check-analytic "${output_file}")
endfunction()

run_checked("${TEST_UTILS}" make-case "${RESULT}" "${CENTERLINE}")
message(STATUS "${COMMAND_OUTPUT}")

set(REFERENCE "${WORK_DIR}/out_omp1.vtp")
run_svslicer(1 "${REFERENCE}")

foreach(repeat RANGE 1 ${REPEATS})
    set(output_file "${WORK_DIR}/out_omp${THREADS}_run${repeat}.vtp")
    run_svslicer(${THREADS} "${output_file}")
    run_checked("${TEST_UTILS}" compare "${REFERENCE}" "${output_file}")
endforeach()

message(STATUS "${REPEATS} runs with ${THREADS} threads match the "
               "single-threaded run and match the analytic solution")

# A runtime thread limit must reduce the number of grid copies as well.
set(ENV{OMP_THREAD_LIMIT} 1)
set(LIMITED "${WORK_DIR}/out_limited.vtp")
run_svslicer(8 "${LIMITED}")
if(NOT SVSLICER_OUTPUT MATCHES "Using 1 thread\\(s\\)")
    message(FATAL_ERROR "Expected one participating thread:\n${SVSLICER_OUTPUT}")
endif()
run_checked("${TEST_UTILS}" compare "${REFERENCE}" "${LIMITED}")
