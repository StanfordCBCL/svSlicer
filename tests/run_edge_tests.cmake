# Check velocity component validation and the thread cap for a single slice.
foreach(var SVSLICER TEST_UTILS WORK_DIR)
    if(NOT DEFINED ${var})
        message(FATAL_ERROR "${var} is not set")
    endif()
endforeach()

file(MAKE_DIRECTORY "${WORK_DIR}")
set(CENTERLINE "${WORK_DIR}/one_point.vtp")
file(WRITE "${CENTERLINE}" [=[
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="1" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">
      <PointData>
        <DataArray type="Float32" Name="CenterlineSectionNormal" NumberOfComponents="3" format="ascii">0 0 1</DataArray>
      </PointData>
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">0.2 0.2 0.25</DataArray></Points>
    </Piece>
  </PolyData>
</VTKFile>
]=])

function(write_case components name tuple result_file)
    string(REPEAT "${tuple} " 4 values)
    file(WRITE "${result_file}" "
<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">
  <UnstructuredGrid>
    <Piece NumberOfPoints=\"4\" NumberOfCells=\"1\">
      <PointData>
        <DataArray type=\"Float64\" Name=\"${name}\" NumberOfComponents=\"${components}\" format=\"ascii\">${values}</DataArray>
        <DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">100 100 100 100</DataArray>
      </PointData>
      <Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">0 0 0  1 0 0  0 1 0  0 0 1</DataArray></Points>
      <Cells>
        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">0 1 2 3</DataArray>
        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">4</DataArray>
        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">10</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>")
endfunction()

foreach(components 1 2)
    set(result_file "${WORK_DIR}/invalid_${components}.vtu")
    if(components EQUAL 1)
        set(tuple "25")
    else()
        set(tuple "0 25")
    endif()
    write_case(${components} VelocityMagnitude "${tuple}" "${result_file}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=8
        "${SVSLICER}" "${result_file}" "${CENTERLINE}" "${WORK_DIR}/invalid_${components}.vtp"
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE output)
    if(NOT result STREQUAL "1" OR
       NOT output MATCHES "VelocityMagnitude.*must have at least 3 components" OR
       output MATCHES "Extract slices")
        message(FATAL_ERROR "Expected velocity validation failure: ${result}\n${output}")
    endif()
endforeach()

function(run_valid result_file output_file threads)
    execute_process(COMMAND ${CMAKE_COMMAND} -E env
        OMP_NUM_THREADS=${threads} OMP_THREAD_LIMIT=8 OMP_DYNAMIC=FALSE
        "${SVSLICER}" "${result_file}" "${CENTERLINE}" "${output_file}"
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE output)
    if(NOT result STREQUAL "0" OR NOT output MATCHES "Using 1 thread\\(s\\)")
        message(FATAL_ERROR "Expected one participating thread: ${result}\n${output}")
    endif()
endfunction()

set(RESULT "${WORK_DIR}/valid.vtu")
set(REFERENCE "${WORK_DIR}/reference.vtp")
write_case(3 Velocity "0 0 25" "${RESULT}")
run_valid("${RESULT}" "${REFERENCE}" 1)
set(OUTPUT "${WORK_DIR}/capped.vtp")
run_valid("${RESULT}" "${OUTPUT}" 8)
execute_process(COMMAND "${TEST_UTILS}" compare "${REFERENCE}" "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
    message(FATAL_ERROR "Capped-thread output differs from the reference")
endif()

write_case(4 Velocity "0 0 25 99" "${RESULT}")
set(OUTPUT "${WORK_DIR}/four_components.vtp")
run_valid("${RESULT}" "${OUTPUT}" 1)
execute_process(COMMAND "${TEST_UTILS}" compare "${REFERENCE}" "${OUTPUT}"
    RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
    message(FATAL_ERROR "Four-component velocity differs from the reference")
endif()
