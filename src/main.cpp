/**
 * @file main.cpp
 * @brief Main svslicer routine.
 *
 */
#include <vtkCellCenters.h>
#include <vtkCellData.h>
#include <vtkCleanPolyData.h>
#include <vtkContourGrid.h>
#include <vtkDoubleArray.h>
#include <vtkExtractCells.h>
#include <vtkKdTree.h>
#include <vtkPointData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <filesystem>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>

/**
 * @brief Check if given string starts with a certain prefix
 *
 * @param string The string to check for the given prefix
 * @param prefix The prefix
 * @return true If string starts with the prefix
 * @return false If string does not start with the prefix
 */
bool starts_with(std::string &string, std::string prefix)
{
    if (string.find(prefix) == 0)
        return true;
    else
        return false;
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
    for (auto const &name : remove_array_names)
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
    for (auto const &name : remove_array_names)
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
    // Create connectivity filter that separates loose parts of the slice
    auto conn_filter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    conn_filter->SetInputData(slice);
    conn_filter->SetExtractionModeToSpecifiedRegions();

    // Minimum distance of a center of a slice region to the position
    double min_distance = 1e6;

    // ID of the current region of the slice filter with connectivity filter
    int region_id = 0;

    // Center of the current slice region
    double center[3] = {0.0, 0.0, 0.0};

    // Clean slice
    vtkSmartPointer<vtkPolyData> clean_slice;
    while (true)
    {
        // Extract next slice region and delete it from the filter
        auto current_slice_region = vtkSmartPointer<vtkPolyData>::New();
        conn_filter->AddSpecifiedRegion(region_id);
        conn_filter->Update();
        current_slice_region->DeepCopy(conn_filter->GetOutput());
        if (current_slice_region->GetNumberOfCells() <= 0)
        {
            break;
        }
        conn_filter->DeleteSpecifiedRegion(region_id);

        // Extract clean poly data
        auto clean_filter = vtkSmartPointer<vtkCleanPolyData>::New();
        clean_filter->SetInputData(current_slice_region);
        clean_filter->Update();
        current_slice_region->DeepCopy(clean_filter->GetOutput());

        // Determine distance between center of all slice points and the
        // position
        auto region_points = current_slice_region->GetPoints();
        int num_comp_points = current_slice_region->GetNumberOfPoints();
        double cx = 0.0;
        double cy = 0.0;
        double cz = 0.0;
        for (int i = 0; i < num_comp_points; i++)
        {
            double pt[3];
            region_points->GetPoint(i, pt);
            cx += pt[0];
            cy += pt[1];
            cz += pt[2];
        }
        center[0] = cx / num_comp_points;
        center[1] = cy / num_comp_points;
        center[2] = cz / num_comp_points;
        double d = (center[0] - position[0]) * (center[0] - position[0]) +
                   (center[1] - position[1]) * (center[1] - position[1]) +
                   (center[2] - position[2]) * (center[2] - position[2]);

        // If current distance is closer than the distances before, keep slice
        // region as clean slice
        if (d < min_distance)
        {
            min_distance = d;
            clean_slice = current_slice_region;
        }
        region_id += 1;
    }
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
                          std::vector<double> &area_of_cells_in_slice)
{
    auto polys = slice->GetPolys();
    double total_area = 0.0;
    int numcells = polys->GetNumberOfCells();
    area_of_cells_in_slice.resize(numcells);

    // Loop over all cells in slice
    for (int icell = 0; icell < numcells; icell++)
    {
        // Extract points of the
        vtkIdType cell_size;
        const vtkIdType *cell_points;
        polys->GetCellAtId(icell, cell_size, cell_points);

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
        slice->GetPoints()->GetPoint(cell_points[0], p0);
        slice->GetPoints()->GetPoint(cell_points[1], p1);
        slice->GetPoints()->GetPoint(cell_points[2], p2);
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
                          std::vector<double> &area_of_cells_in_slice,
                          std::function<double(vtkIdType)> fun)
{
    // Sum of the integrated quantity of interest for each cell
    double sum = 0;

    // Loop over all cells
    for (int j = 0; j < slice->GetPolys()->GetNumberOfCells(); j++)
    {
        // Extract the points belonging to the selected cell
        vtkIdType cell_size;
        const vtkIdType *cell_points;
        slice->GetPolys()->GetCellAtId(j, cell_size, cell_points);

        // Determine mean of quantity of all cell points
        double cur_sum = 0.0;
        for (int k = 0; k < cell_size; k++) cur_sum += fun(cell_points[k]);

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
int main(int argc, char *argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    if (argc < 3)
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
    auto plane_dist = vtkSmartPointer<vtkDoubleArray>::New();
    plane_dist->SetName("plane_dist");
    plane_dist->SetNumberOfComponents(1);
    plane_dist->SetNumberOfTuples(num_pts);
    for (int i = 0; i < num_pts; i++)
    {
        plane_dist->SetValue(i, 0.0);
    }
    result_unstructured_grid->GetPointData()->AddArray(plane_dist);
    result_unstructured_grid->GetPointData()->SetActiveScalars("plane_dist");

    // Get centerline data.
    auto centerline_points = centerline_polydata->GetPoints();
    int num_centerline_points = centerline_polydata->GetNumberOfPoints();
    auto centerline_normals = vtkDoubleArray::SafeDownCast(
        centerline_polydata->GetPointData()->GetArray(
            "CenterlineSectionNormal"));
    auto centerline_max_size = vtkDoubleArray::SafeDownCast(
        centerline_polydata->GetPointData()->GetArray(
            "CenterlineSectionMaxSize"));

    // Prepare area array
    auto area = vtkSmartPointer<vtkDoubleArray>::New();
    area->SetName("area");
    area->SetNumberOfComponents(1);
    area->SetNumberOfTuples(num_centerline_points);
    centerline_polydata->GetPointData()->AddArray(area);

    // Prepare pressure and velocity arrays
    int num_point_arrays =
        result_unstructured_grid->GetPointData()->GetNumberOfArrays();
    for (int ipoint = 0; ipoint < num_point_arrays; ipoint++)
    {
        auto name = std::string(
            result_unstructured_grid->GetPointData()->GetArrayName(ipoint));
        if (starts_with(name, "velocity") || starts_with(name, "pressure"))
        {
            auto array = vtkSmartPointer<vtkDoubleArray>::New();
            array->SetName(name.c_str());
            array->SetNumberOfComponents(1);
            array->SetNumberOfTuples(num_centerline_points);
            centerline_polydata->GetPointData()->AddArray(array);
        }
    }

    // Extract centers of cells and build KdTree to search for cells in radius
    auto cell_center_filter = vtkSmartPointer<vtkCellCenters>::New();
    cell_center_filter->SetInputData(result_unstructured_grid);
    cell_center_filter->Update();
    auto cell_center_kdtree = vtkSmartPointer<vtkKdTree>::New();
    cell_center_kdtree->BuildLocatorFromPoints(
        cell_center_filter->GetOutput()->GetPoints());

    // Start parallel extraction
    std::cout << "Extract slices" << std::endl;
    int num_slices_processed = 0;
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < 10; i++)
    {
        // Extract thread information
        int my_slice;
#pragma omp critical(status)
        {
            num_slices_processed++;
            my_slice = num_slices_processed;
        }

        // Extract current position
        double position1[3];
        double normal1[3];
        double maxsize1[1];
        centerline_points->GetPoint(i, position1);
        centerline_normals->GetTuple(i, normal1);
        centerline_max_size->GetTuple(i, maxsize1);

        // Extract auxilliary second position for interpolation
        double normal2[3];
        double maxsize2[0];
        double position2[3];
        double position2a[3];
        double position2b[3];
        centerline_points->GetPoint(i + 1, position2a);
        centerline_points->GetPoint(i - 1, position2b);

        double quaddista = 0.0;
        double quaddistb = 0.0;
        for (int j = 0; j < 3; j++)
        {
            quaddista +=
                (position1[j] - position2a[j]) * (position1[j] - position2a[j]);
            quaddistb +=
                (position1[j] - position2b[j]) * (position1[j] - position2b[j]);
        }

        if (quaddista < quaddistb)
        {
            position2[0] = position2a[0];
            position2[1] = position2a[1];
            position2[2] = position2a[2];
            centerline_normals->GetTuple(i + 1, normal2);
            centerline_max_size->GetTuple(i + 1, maxsize2);
        }
        else
        {
            position2[0] = position2b[0];
            position2[1] = position2b[1];
            position2[2] = position2b[2];
            centerline_normals->GetTuple(i - 1, normal2);
            centerline_max_size->GetTuple(i - 1, maxsize2);
        }

        // Setup position and slice
        vtkSmartPointer<vtkPolyData> slice;
        vtkSmartPointer<vtkUnstructuredGrid> my_grid;
        double position[3];
        double normal[3];
        double srad;

        // Start slice extraction routine
        float weight = 0.01;
        while (true)
        {
            // Determine interpolated position, normal and max_size
            double nnorm = 0.0;
            for (int j = 0; j < 3; j++)
            {
                position[j] =
                    position1[j] * (1.0 - weight) + position2[j] * weight;
                normal[j] = normal1[j] * (1.0 - weight) + normal2[j] * weight;
                nnorm += normal[j] * normal[j];
            }
            for (int j = 0; j < 3; j++)
            {
                normal[j] = normal[j] / nnorm;
            }
            srad = (maxsize1[0] * (1.0 - weight) + maxsize2[0] * weight) * 0.55;

            // Find cells of interest by searching in radius around position
            auto selected_cells = vtkSmartPointer<vtkIdList>::New();
            cell_center_kdtree->FindPointsWithinRadius(srad, position,
                                                       selected_cells);
            auto cell_extractor = vtkSmartPointer<vtkExtractCells>::New();
            cell_extractor->SetInputData(result_unstructured_grid);
            cell_extractor->SetCellList(selected_cells);
            cell_extractor->Update();
            my_grid = cell_extractor->GetOutput();

            // Compute distance of each mesh point from the slicing plane.
            for (int iter1 = 0; iter1 < my_grid->GetNumberOfPoints(); iter1++)
            {
                double pt[3];
                my_grid->GetPoints()->GetPoint(iter1, pt);
                double d[1];
                d[0] = normal[0] * (pt[0] - position[0]) +
                       normal[1] * (pt[1] - position[1]) +
                       normal[2] * (pt[2] - position[2]);
                my_grid->GetPointData()
                    ->GetArray("plane_dist")
                    ->SetTuple(iter1, d);
            }

            // Set the scalar field used to extract a slice.
            my_grid->GetPointData()->SetActiveScalars("plane_dist");

            // Extract slice with contour filter where distance to plane is zero
            auto contour = vtkSmartPointer<vtkContourGrid>::New();
            contour->SetInputData(my_grid);
            contour->SetValue(0, 0.0);
            contour->ComputeNormalsOff();
            contour->Update();
            slice = contour->GetOutput();

            // Break if extracted slice is not emtpy
            if (slice->GetNumberOfPoints())
            {
                break;
            }

            // Raise runtime error if no slice found after 5 iterations
            if (weight == 0.6)
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
        int numcells = slice->GetPolys()->GetNumberOfCells();

        // Calculate area of slice and of all cells in slice
        std::vector<double> area_of_cells_in_slice;
        double area_i[1];
        area_i[0] = compute_slice_area(slice, area_of_cells_in_slice);
        centerline_polydata->GetPointData()->GetArray("area")->SetTuple(i,
                                                                        area_i);
        double area_inv = 1.0 / area_i[0];

        // Iterate over point data arrays and integrate over slice
        for (int ipoint = 0; ipoint < num_point_arrays_slice; ipoint++)
        {
            auto name =
                std::string(slice->GetPointData()->GetArrayName(ipoint));

            // Integrate velocity profile
            if (starts_with(name, "velocity"))
            {
                // Calculate flux (flow in normal direction)
                double flux[num_points_slice];
                for (int j = 0; j < num_points_slice; j++)
                {
                    double *vel =
                        slice->GetPointData()->GetArray(ipoint)->GetTuple3(j);
                    double curflux = 0;
                    for (int k = 0; k < 3; k++)
                    {
                        curflux += vel[k] * normal[k];
                    }
                    flux[j] = curflux;
                }

                // Integrate flux over slice
                double result[1];
                result[0] = integrate_on_slice(slice, area_of_cells_in_slice,
                                               [&](vtkIdType index) -> double
                                               { return flux[index]; });
                centerline_polydata->GetPointData()
                    ->GetArray(name.c_str())
                    ->SetTuple(i, result);
            }

            // Integrate pressure
            if (starts_with(name, "pressure"))
            {
                // Integrate pressure over slice area
                double result[1];
                result[0] = integrate_on_slice(slice, area_of_cells_in_slice,
                                               [&](vtkIdType index) -> double {
                                                   return slice->GetPointData()
                                                       ->GetArray(ipoint)
                                                       ->GetTuple1(index);
                                               });
                result[0] *= area_inv;
                centerline_polydata->GetPointData()
                    ->GetArray(name.c_str())
                    ->SetTuple(i, result);
            }
        }
        std::cout << "Completed slice " << my_slice << "/"
                  << num_centerline_points << std::endl;
    }

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = stop - start;
    std::cout << "Slice extraction completed in " << duration.count() << " s ("
              << duration.count() / num_slices_processed << " s per slice)"
              << std::endl;
}
