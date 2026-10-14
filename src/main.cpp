/**
 * @file main.cpp
 * @brief Main svslicer routine.
 *
 */
#include <omp.h>
#include <vtkCellData.h>
#include <vtkCleanPolyData.h>
#include <vtkCutter.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkSMPTools.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersionMacros.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @brief Check if given string starts with a certain prefix, ignoring case
 *
 * Case is ignored so that both svSolver ("velocity", "pressure") and
 * svMultiPhysics ("Velocity", "Pressure") result arrays are recognized.
 *
 * @param string The string to check for the given prefix
 * @param prefix The prefix
 * @return true If string starts with the prefix
 * @return false If string does not start with the prefix
 */
bool starts_with(const std::string& string, const std::string& prefix)
{
    if (string.size() < prefix.size())
    {
        return false;
    }
    for (size_t i = 0; i < prefix.size(); i++)
    {
        if (std::tolower(static_cast<unsigned char>(string[i])) !=
            std::tolower(static_cast<unsigned char>(prefix[i])))
        {
            return false;
        }
    }
    return true;
}

/**
 * @brief Clean-up point data arrays.
 *
 * Delete all point data arrays that are not marked for keeping.
 *
 * @tparam T Type of data (e.g. vtkPolyData, vtkUnstructuredGrid)
 * @param data_pointer VTK smart pointer to the data
 * @param keep_names Set of array names to keep
 */
template <typename T>
void clean_up_point_data(vtkSmartPointer<T> data_pointer,
                         std::set<std::string> keep_names)
{
    // Get number of point data arrays
    int num_point_arrays = data_pointer->GetPointData()->GetNumberOfArrays();

    // Vector storing the names of the arrays to remove
    std::vector<std::string> remove_array_names;

    // Iterate over all point data arrays
    for (int i = 0; i < num_point_arrays; i++)
    {
        auto array_name =
            std::string(data_pointer->GetPointData()->GetArrayName(i));

        // Mark array for removal if name not in keep_names
        bool keep = false;
        for (auto name : keep_names)
        {
            if (starts_with(array_name, name))
            {
                keep = true;
            }
        }

        // Add array to the list of data arrays for removal if it should not be
        // kept
        if (keep)
        {
            std::cout << "Keeping point data array '" << array_name << "'"
                      << std::endl;
        }
        else
        {
            std::cout << "Removing point data array '" << array_name << "'"
                      << std::endl;
            remove_array_names.push_back(array_name);
        }
    }

    // Delete all data arrays marked for removal
    for (auto const& name : remove_array_names)
    {
        data_pointer->GetPointData()->RemoveArray(name.c_str());
    }
}

/**
 * @brief Clean-up cell data arrays.
 *
 * Delete all cell data arrays that are not marked for keeping.
 *
 * @tparam T Type of data (e.g. vtkPolyData, vtkUnstructuredGrid)
 * @param data_pointer VTK smart pointer to the data
 * @param keep_names Set of array names to keep
 */
template <typename T>
void clean_up_cell_data(vtkSmartPointer<T> data_pointer,
                        std::set<std::string> keep_names)
{
    // Get number of cell data arrays
    int num_cell_arrays = data_pointer->GetCellData()->GetNumberOfArrays();

    // Vector storing the names of the arrays to remove
    std::vector<std::string> remove_array_names;

    // Iterate over all cell data arrays
    for (int i = 0; i < num_cell_arrays; i++)
    {
        auto array_name =
            std::string(data_pointer->GetCellData()->GetArrayName(i));

        // Mark array for removal if name not in keep_names
        bool keep = false;
        for (auto name : keep_names)
        {
            if (starts_with(array_name, name))
            {
                keep = true;
            }
        }

        // Add array to the list of data arrays for removal if it should not be
        // kept
        if (keep)
        {
            std::cout << "Keeping cell data array '" << array_name << "'"
                      << std::endl;
        }
        else
        {
            std::cout << "Removing cell data array '" << array_name << "'"
                      << std::endl;
            remove_array_names.push_back(array_name);
        }
    }

    // Delete all data arrays marked for removal
    for (auto const& name : remove_array_names)
    {
        data_pointer->GetCellData()->RemoveArray(name.c_str());
    }
}

/**
 * @brief Clean up the slice to only contain region of interest
 *
 * Extracts the slice region of interest. In case a slice has cut through a
 * volume multiple times, this function extracts only the part of the slice
 * that is of interest.
 *
 * @param position The centerline position, which the slice originated from
 * @param slice The slice
 * @return vtkSmartPointer<vtkPolyData> The new slice
 */
vtkSmartPointer<vtkPolyData> get_clean_slice(double position[3],
                                             vtkSmartPointer<vtkPolyData> slice)
{
    // Extract clean poly data
    auto clean_filter = vtkSmartPointer<vtkCleanPolyData>::New();
    clean_filter->SetInputData(slice);
    clean_filter->Update();

    // Create connectivity filter that separates loose parts of the slice
    auto conn_filter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    conn_filter->SetInputData(clean_filter->GetOutput());
    conn_filter->SetExtractionModeToClosestPointRegion();
    conn_filter->SetClosestPoint(position);
    auto clean_slice = vtkSmartPointer<vtkPolyData>::New();
    conn_filter->Update();
    clean_slice->DeepCopy(conn_filter->GetOutput());

    return clean_slice;
}

/**
 * @brief Compute the area of a given slice and its individual cells
 *
 * @param slice The slice to compute the area of
 * @param area_of_cells_in_slice The array to save the area of the cells in
 * @return double The total area of the slice.
 */
double compute_slice_area(vtkSmartPointer<vtkPolyData> slice,
                          std::vector<double>& area_of_cells_in_slice)
{
    auto polys = slice->GetPolys();
    double total_area = 0.0;
    int numcells = polys->GetNumberOfCells();
    area_of_cells_in_slice.resize(numcells);
    auto cell_points = vtkSmartPointer<vtkIdList>::New();

    // Loop over all cells in slice. The sequential traversal is used instead
    // of GetCellAtId, which is not available in VTK < 9.
    polys->InitTraversal();
    for (int icell = 0; icell < numcells; icell++)
    {
        // Extract points of the cell
        if (!polys->GetNextCell(cell_points))
        {
            throw std::runtime_error("Failed to traverse slice cells.");
        }
        vtkIdType cell_size = cell_points->GetNumberOfIds();

        // From now on, we assume that the elements are triangles (true for
        // tetrahedral mesh)
        if (cell_size != 3)
        {
            throw std::runtime_error("Only triangles supported in slice.");
        }

        // Calculate area of triangle cell
        double p0[3];
        double p1[3];
        double p2[3];
        slice->GetPoints()->GetPoint(cell_points->GetId(0), p0);
        slice->GetPoints()->GetPoint(cell_points->GetId(1), p1);
        slice->GetPoints()->GetPoint(cell_points->GetId(2), p2);
        double cell_area = vtkTriangle::TriangleArea(p0, p1, p2);

        // Add cell area to list of cell areas and add to total slice area
        area_of_cells_in_slice[icell] = cell_area;
        total_area += cell_area;
    }
    return total_area;
}

/**
 * @brief Integrate point data quantity over slice
 *
 * @param slice The slice to integrate quantity
 * @param area_of_cells_in_slice The areas of the cells in the slice
 * @param fun The function to extract the quantity for a given point ID
 * @return double The integrated quantity
 */
double integrate_on_slice(vtkSmartPointer<vtkPolyData> slice,
                          std::vector<double>& area_of_cells_in_slice,
                          std::function<double(vtkIdType)> fun)
{
    // Sum of the integrated quantity of interest for each cell
    double sum = 0;
    auto polys = slice->GetPolys();
    auto cell_points = vtkSmartPointer<vtkIdList>::New();

    // Loop over all cells
    polys->InitTraversal();
    for (int j = 0; j < polys->GetNumberOfCells(); j++)
    {
        // Extract the points belonging to the selected cell
        if (!polys->GetNextCell(cell_points))
        {
            throw std::runtime_error("Failed to traverse slice cells.");
        }
        vtkIdType cell_size = cell_points->GetNumberOfIds();

        // Determine mean of quantity of all cell points
        double cur_sum = 0.0;
        for (vtkIdType k = 0; k < cell_size; k++)
            cur_sum += fun(cell_points->GetId(k));

        // Add cell value times the cell area to the total sum
        sum += area_of_cells_in_slice[j] * (cur_sum) / cell_size;
    }
    return sum;
}

/**
 * @brief Main routine of svSlicer.
 *
 * Performs the slicing of 3D volumetric results to map the quantities of
 * interest (pressure and flow/velocity) to a given centerline.
 *
 * @param argc Number of command line arguments
 * @param argv Command line arguments passed to svSlicer
 * @return int
 */
int main(int argc, char* argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    if (argc < 4)
    {
        std::cerr
            << "Usage: svslicer [RESULT_FILE] [CENTERLINE_FILE] [OUTPUT_FILE]"
            << std::endl;
        return EXIT_FAILURE;
    }

    // Check if result file exists
    std::string result_file_name = argv[1];
    if (std::filesystem::exists(result_file_name) == false)
    {
        std::cerr << "ERROR: 3D result file '" << result_file_name
                  << "' does not exist" << std::endl;
        return EXIT_FAILURE;
    }

    // Check if centerline file exists
    std::string centerline_file_name = argv[2];
    if (std::filesystem::exists(centerline_file_name) == false)
    {
        std::cerr << "ERROR: Centerline file '" << centerline_file_name
                  << "' does not exist" << std::endl;
        return EXIT_FAILURE;
    }

    // Read in the result file.
    std::cout << "Reading result file '" << result_file_name << "'"
              << std::endl;
    auto result_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    result_reader->SetFileName(result_file_name.c_str());
    result_reader->Update();
    auto result_unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    result_unstructured_grid = result_reader->GetOutput();
    clean_up_point_data<vtkUnstructuredGrid>(result_unstructured_grid,
                                             {"pressure", "velocity"});
    clean_up_cell_data<vtkUnstructuredGrid>(result_unstructured_grid,
                                            {"pressure", "velocity"});

    // Read in the centerline file.
    std::cout << "Reading centerline file '" << centerline_file_name << "'"
              << std::endl;
    auto centerline_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    centerline_reader->SetFileName(centerline_file_name.c_str());
    centerline_reader->Update();
    auto centerline_polydata = vtkSmartPointer<vtkPolyData>::New();
    centerline_polydata = centerline_reader->GetOutput();
    clean_up_point_data<vtkPolyData>(
        centerline_polydata,
        {"BifurcationId", "BifurcationIdTmp", "BranchId", "BranchIdTmp",
         "CenterlineId", "CenterlineSectionArea",
         "CenterlineSectionBifurcation", "CenterlineSectionClosed",
         "CenterlineSectionMaxSize", "CenterlineSectionMinSize",
         "CenterlineSectionNormal", "CenterlineSectionRemove",
         "CenterlineSectionShape", "GlobalNodeId",
         "MaximumInscribedSphereRadius", "Path"});
    clean_up_cell_data<vtkPolyData>(centerline_polydata, {});

    // Add a point data array to store plane distance.
    int num_pts = result_unstructured_grid->GetNumberOfPoints();

    // Get centerline data.
    int num_centerline_points = centerline_polydata->GetNumberOfPoints();
    auto centerline_points = centerline_polydata->GetPoints();
    // The normals are read through the generic vtkDataArray interface, so
    // that both Float32 and Float64 arrays are accepted.
    vtkDataArray* centerline_normals =
        centerline_polydata->GetPointData()->GetArray(
            "CenterlineSectionNormal");
    if (centerline_normals == nullptr ||
        centerline_normals->GetNumberOfComponents() != 3)
    {
        std::cerr << "ERROR: Centerline file '" << centerline_file_name
                  << "' has no 3-component 'CenterlineSectionNormal' point "
                     "data array"
                  << std::endl;
        return EXIT_FAILURE;
    }

    // Prepare pressure and velocity arrays. The threads write their results
    // directly into the raw array memory (each thread a different index)
    // rather than through VTK accessors.
    std::map<std::string, double*> output_values;
    int num_point_arrays =
        result_unstructured_grid->GetPointData()->GetNumberOfArrays();
    for (int point_array_idx = 0; point_array_idx < num_point_arrays;
         point_array_idx++)
    {
        auto name =
            std::string(result_unstructured_grid->GetPointData()->GetArrayName(
                point_array_idx));
        if (starts_with(name, "velocity") &&
            result_unstructured_grid->GetPointData()
                    ->GetArray(point_array_idx)
                    ->GetNumberOfComponents() < 3)
        {
            std::cerr << "ERROR: Velocity array '" << name << "' in '"
                      << result_file_name << "' must have at least 3 components"
                      << std::endl;
            return EXIT_FAILURE;
        }
        if (starts_with(name, "velocity") || starts_with(name, "pressure"))
        {
            auto array = vtkSmartPointer<vtkDoubleArray>::New();
            array->SetName(name.c_str());
            array->SetNumberOfComponents(1);
            array->SetNumberOfTuples(num_centerline_points);
            centerline_polydata->GetPointData()->AddArray(array);
            output_values[name] = array->GetPointer(0);
        }
    }

    // vtkCutter must not be run concurrently on a shared input. VTK 7.1
    // evaluates the cut function through vtkDataSet::GetPoint(vtkIdType),
    // which returns a pointer into a tuple buffer that belongs to the points
    // array and is shared by all callers, and every cutter also writes
    // pipeline information and cached bounds into its input. Give every
    // thread its own deep copy of the grid (thread 0 uses the original), so
    // that no VTK object is shared between threads. A shallow copy is not
    // enough, because it shares the points and point data arrays.
    int num_threads =
        std::max(1, std::min(omp_get_max_threads(), num_centerline_points));

    // The slices are processed in parallel by OpenMP, so VTK must not run
    // its own SMP threads (TBB or std::thread in some VTK 9 builds) inside
    // each cut. They oversubscribe the cores, and their scheduling changes
    // the order of the slice triangles and with it the floating-point sums.
    // Limiting TBB to one thread is not enough, because other OpenMP threads
    // can still execute parts of a cut, so use the Sequential backend where
    // it can be selected at run time (VTK >= 9.1).
#if VTK_MAJOR_VERSION > 9 || (VTK_MAJOR_VERSION == 9 && VTK_MINOR_VERSION >= 1)
    if (!vtkSMPTools::SetBackend("Sequential"))
    {
        std::cerr << "WARNING: could not select VTK's Sequential SMP backend; "
                  << "using " << vtkSMPTools::GetBackend() << std::endl;
    }
#endif
    vtkSMPTools::Initialize(1);
    std::vector<vtkSmartPointer<vtkUnstructuredGrid>> thread_grids;

    // Start parallel extraction
    std::cout << "Extract slices" << std::endl;
    auto slice_start = std::chrono::high_resolution_clock::now();
    int num_slices_processed = 0;
    // Exceptions must not leave the OpenMP region, which would call
    // std::terminate. The first error is recorded, the remaining slices are
    // skipped, and the error is reported after the loop.
    std::atomic<bool> slice_failed(false);
    std::string slice_error;
#pragma omp parallel num_threads(num_threads)
    {
#pragma omp single
        {
            // The runtime may reduce the requested team size. Copy only for
            // the actual team, serially, before any thread starts slicing.
            try
            {
                int team_size = omp_get_num_threads();
                std::cout << "Using " << team_size << " thread(s), "
                          << result_unstructured_grid->GetActualMemorySize() /
                                 1024.0
                          << " MB per copy of the 3D result" << std::endl;
                thread_grids.resize(team_size);
                thread_grids[0] = result_unstructured_grid;
                for (int t = 1; t < team_size; t++)
                {
                    thread_grids[t] =
                        vtkSmartPointer<vtkUnstructuredGrid>::New();
                    thread_grids[t]->DeepCopy(result_unstructured_grid);
                }
            }
            catch (const std::exception& e)
            {
                slice_error =
                    std::string("Failed to prepare thread grids: ") + e.what();
                slice_failed = true;
            }
            catch (...)
            {
                slice_error = "Failed to prepare thread grids: unknown error";
                slice_failed = true;
            }
        }
#pragma omp for schedule(dynamic)
        for (int i = 0; i < num_centerline_points; i++)
        {
            if (slice_failed)
            {
                continue;
            }
            try
            {
                // Extract thread information
                int my_slice;
#pragma omp critical(status)
                {
                    num_slices_processed++;
                    my_slice = num_slices_processed;
                }

                // Grid owned by this thread
                vtkUnstructuredGrid* grid = thread_grids[omp_get_thread_num()];

                // Extract current position. The centerline is shared, but only
                // read through the accessors that copy into a caller-provided
                // buffer.
                double position1[3];
                double normal1[3];
                centerline_points->GetPoint(i, position1);
                centerline_normals->GetTuple(i, normal1);

                // Extract auxilliary second position for interpolation. The
                // neighbor that is closest to the current point is used. Guard
                // against the first/last centerline point, which only has a
                // single neighbor.
                double normal2[3];
                double position2[3];
                double position2a[3];
                double position2b[3];
                bool has_next = (i + 1) < num_centerline_points;
                bool has_prev = (i - 1) >= 0;

                double quaddista = std::numeric_limits<double>::max();
                double quaddistb = std::numeric_limits<double>::max();
                if (has_next)
                {
                    centerline_points->GetPoint(i + 1, position2a);
                    quaddista = 0.0;
                    for (int j = 0; j < 3; j++)
                        quaddista += (position1[j] - position2a[j]) *
                                     (position1[j] - position2a[j]);
                }
                if (has_prev)
                {
                    centerline_points->GetPoint(i - 1, position2b);
                    quaddistb = 0.0;
                    for (int j = 0; j < 3; j++)
                        quaddistb += (position1[j] - position2b[j]) *
                                     (position1[j] - position2b[j]);
                }

                if (has_next && quaddista <= quaddistb)
                {
                    position2[0] = position2a[0];
                    position2[1] = position2a[1];
                    position2[2] = position2a[2];
                    centerline_normals->GetTuple(i + 1, normal2);
                }
                else if (has_prev)
                {
                    position2[0] = position2b[0];
                    position2[1] = position2b[1];
                    position2[2] = position2b[2];
                    centerline_normals->GetTuple(i - 1, normal2);
                }
                else
                {
                    // Degenerate single-point centerline: no neighbor to
                    // interpolate towards, so fall back to the point itself.
                    position2[0] = position1[0];
                    position2[1] = position1[1];
                    position2[2] = position1[2];
                    normal2[0] = normal1[0];
                    normal2[1] = normal1[1];
                    normal2[2] = normal1[2];
                }

                // Setup slice
                vtkSmartPointer<vtkPolyData> slice;
                double position[3];
                double normal[3];

                // Start slice extraction routine
                float weight = 0.0;
                while (true)
                {
                    // Determine interpolated position, normal and max_size
                    double nnorm = 0.0;
                    for (int j = 0; j < 3; j++)
                    {
                        position[j] = position1[j] * (1.0 - weight) +
                                      position2[j] * weight;
                        normal[j] =
                            normal1[j] * (1.0 - weight) + normal2[j] * weight;
                        nnorm += normal[j] * normal[j];
                    }
                    nnorm = std::sqrt(nnorm);
                    for (int j = 0; j < 3; j++)
                    {
                        normal[j] = normal[j] / nnorm;
                    }

                    auto plane = vtkSmartPointer<vtkPlane>::New();
                    plane->SetOrigin(position);
                    plane->SetNormal(normal);

                    auto cutter = vtkSmartPointer<vtkCutter>::New();
                    cutter->SetCutFunction(plane);
                    cutter->SetInputData(grid);
                    cutter->Update();
                    slice = cutter->GetOutput();

                    // Update weight
                    weight += 0.1;

                    // Break if extracted slice is not emtpy
                    if (slice->GetNumberOfPoints())
                    {
                        break;
                    }
                    std::cout << "Slice was emtpy." << std::endl;

                    // Raise runtime error if no slice found after 5 iterations
                    if (weight >= 0.6)
                    {
                        throw std::runtime_error("Interpolation failed.");
                    }
                }

                // Clean up the slice
                slice = get_clean_slice(position, slice);

                // Extract slice data
                vtkIdType num_point_arrays_slice =
                    slice->GetPointData()->GetNumberOfArrays();
                int num_points_slice = slice->GetNumberOfPoints();

                // Calculate area of slice and of all cells in slice
                std::vector<double> area_of_cells_in_slice;
                double slice_area;
                slice_area = compute_slice_area(slice, area_of_cells_in_slice);
                double slice_area_inv = 1.0 / slice_area;

                // Iterate over point data arrays and integrate over slice
                for (int point_array_idx = 0;
                     point_array_idx < num_point_arrays_slice;
                     point_array_idx++)
                {
                    auto name = std::string(
                        slice->GetPointData()->GetArrayName(point_array_idx));

                    // Integrate velocity profile
                    if (starts_with(name, "velocity"))
                    {
                        // Calculate flux (flow in normal direction)
                        std::vector<double> flux(num_points_slice);
                        for (int j = 0; j < num_points_slice; j++)
                        {
                            // Read the first three components one at a time, so
                            // an array with more components cannot overflow a
                            // buffer.
                            vtkDataArray* vel = slice->GetPointData()->GetArray(
                                point_array_idx);
                            double curflux = 0.0;
                            for (int k = 0; k < 3; k++)
                            {
                                curflux += vel->GetComponent(j, k) * normal[k];
                            }
                            flux[j] = curflux;
                        }

                        // Integrate flux over slice
                        output_values.at(name)[i] =
                            integrate_on_slice(slice, area_of_cells_in_slice,
                                               [&](vtkIdType index) -> double
                                               { return flux[index]; });
                    }

                    // Integrate pressure
                    if (starts_with(name, "pressure"))
                    {
                        // Integrate pressure over slice area
                        double result = integrate_on_slice(
                            slice, area_of_cells_in_slice,
                            [&](vtkIdType index) -> double
                            {
                                return slice->GetPointData()
                                    ->GetArray(point_array_idx)
                                    ->GetTuple1(index);
                            });
                        output_values.at(name)[i] = result * slice_area_inv;
                    }
                }
                std::cout << "Completed slice " << my_slice << "/"
                          << num_centerline_points << std::endl;
            }
            catch (const std::exception& e)
            {
#pragma omp critical(error)
                {
                    if (!slice_failed)
                    {
                        slice_error = "Slice at centerline point " +
                                      std::to_string(i) +
                                      " failed: " + e.what();
                        slice_failed = true;
                    }
                }
            }
            catch (...)
            {
#pragma omp critical(error)
                {
                    if (!slice_failed)
                    {
                        slice_error = "Slice at centerline point " +
                                      std::to_string(i) +
                                      " failed with an unknown error";
                        slice_failed = true;
                    }
                }
            }
        }
    }
    if (slice_failed)
    {
        std::cerr << "ERROR: " << slice_error << std::endl;
        return EXIT_FAILURE;
    }
    auto slice_stop = std::chrono::high_resolution_clock::now();

    std::string output_file_name = argv[3];
    std::cout << "Writing output file '" << output_file_name << "'"
              << std::endl;
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(output_file_name.c_str());
    writer->SetInputData(centerline_polydata);
    writer->Update();
    writer->Write();

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = stop - start;
    std::chrono::duration<double> duration_sice = slice_stop - slice_start;
    std::cout << "Slice extraction completed in " << duration.count() << " s ("
              << duration_sice.count() / num_centerline_points
              << " s per slice)" << std::endl;
}
