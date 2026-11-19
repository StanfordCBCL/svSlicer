/**
 * @file svslicer_test_utils.cpp
 * @brief Test helpers for svSlicer: synthetic test case and output checks.
 *
 * Usage:
 *   svslicer_test_utils make-case RESULT.vtu CENTERLINE.vtp
 *   svslicer_test_utils check-analytic OUTPUT.vtp
 *   svslicer_test_utils compare REFERENCE.vtp OUTPUT.vtp
 *
 * The test case is a straight tube along z, meshed with linear tetrahedra
 * (a staircase approximation of a circular cross-section, so that every
 * slice z = const has an exactly known area). The velocity and pressure
 * fields are linear, so linear interpolation on the slice reproduces them
 * exactly and the expected flow and mean pressure are known analytically.
 */
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSMPTools.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <string>

namespace
{
// Tube geometry and mesh resolution
const double radius = 1.0;
const double length = 4.0;
const int num_cells_xy = 40;
const int num_cells_z = 80;
const double h_xy = 2.0 * radius / num_cells_xy;
const double h_z = length / num_cells_z;

// Linear fields: velocity_z = U (1 + a x + b y), pressure = p0 - g z
const double mean_velocity = 25.0;
const double vel_x_slope = 0.3;
const double vel_y_slope = 0.2;
const double pressure_inlet = 13332.2;
const double pressure_gradient = 400.0;

// Relative tolerance for the analytic check
const double analytic_tolerance = 1e-9;

// Relative tolerance between the multithreaded and single-threaded results.
// It allows only summation-order rounding (a few 1e-16); the race on a
// shared input grid produced errors of 1e-3 and more.
const double compare_tolerance = 1e-12;

/**
 * @brief Check if the column of hexahedra at (i, j) is inside the tube
 */
bool column_inside(int i, int j)
{
    double x = -radius + (i + 0.5) * h_xy;
    double y = -radius + (j + 0.5) * h_xy;
    return x * x + y * y < radius * radius;
}

/**
 * @brief Cross-sectional area of the staircase tube
 */
double tube_area()
{
    int num_columns = 0;
    for (int i = 0; i < num_cells_xy; i++)
        for (int j = 0; j < num_cells_xy; j++)
            if (column_inside(i, j)) num_columns++;
    return num_columns * h_xy * h_xy;
}

vtkIdType lattice_id(int i, int j, int k)
{
    return (static_cast<vtkIdType>(k) * (num_cells_xy + 1) + j) *
               (num_cells_xy + 1) +
           i;
}

int make_case(const std::string& result_file,
              const std::string& centerline_file)
{
    // Lattice points and linear fields
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToDouble();
    auto velocity = vtkSmartPointer<vtkDoubleArray>::New();
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(3);
    auto pressure = vtkSmartPointer<vtkDoubleArray>::New();
    pressure->SetName("Pressure");
    auto unused = vtkSmartPointer<vtkDoubleArray>::New();
    unused->SetName("WSS");
    for (int k = 0; k <= num_cells_z; k++)
        for (int j = 0; j <= num_cells_xy; j++)
            for (int i = 0; i <= num_cells_xy; i++)
            {
                double x = -radius + i * h_xy;
                double y = -radius + j * h_xy;
                double z = k * h_z;
                points->InsertNextPoint(x, y, z);
                velocity->InsertNextTuple3(
                    0.1 * y, -0.1 * x,
                    mean_velocity * (1.0 + vel_x_slope * x + vel_y_slope * y));
                pressure->InsertNextValue(pressure_inlet -
                                          pressure_gradient * z);
                unused->InsertNextValue(1.0);
            }

    // Split every hexahedron inside the tube into six tetrahedra along its
    // main diagonal (Kuhn triangulation, conforming between neighbors).
    const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                             {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    for (int k = 0; k < num_cells_z; k++)
        for (int j = 0; j < num_cells_xy; j++)
            for (int i = 0; i < num_cells_xy; i++)
            {
                if (!column_inside(i, j)) continue;
                for (auto const& perm : perms)
                {
                    int corner[3] = {0, 0, 0};
                    vtkIdType ids[4];
                    ids[0] = lattice_id(i, j, k);
                    for (int v = 0; v < 3; v++)
                    {
                        corner[perm[v]] = 1;
                        ids[v + 1] = lattice_id(i + corner[0], j + corner[1],
                                                k + corner[2]);
                    }
                    cells->InsertNextCell(4, ids);
                }
            }

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    grid->SetCells(VTK_TETRA, cells);
    grid->GetPointData()->AddArray(velocity);
    grid->GetPointData()->AddArray(pressure);
    grid->GetPointData()->AddArray(unused);

    auto grid_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    grid_writer->SetFileName(result_file.c_str());
    // Base64-encode the data: the VTK 9.1 XML reader shipped with Ubuntu
    // 24.04 fails on large raw appended data with expat >= 2.6.
    grid_writer->SetDataModeToBinary();
    grid_writer->SetInputData(grid);
    if (!grid_writer->Write()) return EXIT_FAILURE;

    // Centerline along the tube axis. The planes are placed between mesh
    // layers so that no slice passes exactly through mesh points. The normals
    // are stored as Float32, as written by many Python tools.
    auto cl_points = vtkSmartPointer<vtkPoints>::New();
    auto normals = vtkSmartPointer<vtkFloatArray>::New();
    normals->SetName("CenterlineSectionNormal");
    normals->SetNumberOfComponents(3);
    auto branch = vtkSmartPointer<vtkDoubleArray>::New();
    branch->SetName("BranchId");
    auto path = vtkSmartPointer<vtkDoubleArray>::New();
    path->SetName("Path");
    auto line = vtkSmartPointer<vtkIdList>::New();
    for (int k = 2; k < num_cells_z - 2; k++)
    {
        double z = (k + 0.37) * h_z;
        line->InsertNextId(cl_points->InsertNextPoint(0.0, 0.0, z));
        normals->InsertNextTuple3(0.0, 0.0, 1.0);
        branch->InsertNextValue(0.0);
        path->InsertNextValue(z);
    }
    auto lines = vtkSmartPointer<vtkCellArray>::New();
    lines->InsertNextCell(line);
    auto centerline = vtkSmartPointer<vtkPolyData>::New();
    centerline->SetPoints(cl_points);
    centerline->SetLines(lines);
    centerline->GetPointData()->AddArray(normals);
    centerline->GetPointData()->AddArray(branch);
    centerline->GetPointData()->AddArray(path);

    auto cl_writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    cl_writer->SetFileName(centerline_file.c_str());
    cl_writer->SetDataModeToBinary();
    cl_writer->SetInputData(centerline);
    if (!cl_writer->Write()) return EXIT_FAILURE;

    std::cout << "VTK " << vtkVersion::GetVTKVersion()
#if VTK_MAJOR_VERSION >= 9
              << ", SMP backend " << vtkSMPTools::GetBackend()
#endif
              << ", " << vtkSMPTools::GetEstimatedNumberOfThreads()
              << " SMP thread(s)" << std::endl;
    std::cout << "Wrote " << grid->GetNumberOfCells() << " tetrahedra and "
              << centerline->GetNumberOfPoints() << " centerline points"
              << std::endl;
    return EXIT_SUCCESS;
}

vtkSmartPointer<vtkPolyData> read_output(const std::string& file_name)
{
    auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(file_name.c_str());
    reader->Update();
    auto output = vtkSmartPointer<vtkPolyData>::New();
    output->DeepCopy(reader->GetOutput());
    for (auto name : {"Velocity", "Pressure"})
    {
        if (output->GetPointData()->GetArray(name) == nullptr)
        {
            std::cerr << "ERROR: '" << file_name << "' has no '" << name
                      << "' array" << std::endl;
            return nullptr;
        }
    }
    return output;
}

int check_analytic(const std::string& output_file)
{
    auto output = read_output(output_file);
    if (!output) return EXIT_FAILURE;
    auto flow = output->GetPointData()->GetArray("Velocity");
    auto pres = output->GetPointData()->GetArray("Pressure");

    // The x and y slopes of the velocity integrate to zero over the
    // symmetric cross-section.
    double expected_flow = mean_velocity * tube_area();
    int num_failed = 0;
    for (vtkIdType i = 0; i < output->GetNumberOfPoints(); i++)
    {
        double x[3];
        output->GetPoint(i, x);
        double expected_pressure = pressure_inlet - pressure_gradient * x[2];
        double flow_err =
            std::abs(flow->GetTuple1(i) - expected_flow) / expected_flow;
        double pres_err = std::abs(pres->GetTuple1(i) - expected_pressure) /
                          expected_pressure;
        if (!(flow_err <= analytic_tolerance) ||
            !(pres_err <= analytic_tolerance))
        {
            if (num_failed < 10)
            {
                std::cerr << "Point " << i << " (z = " << x[2] << "): flow "
                          << flow->GetTuple1(i) << " (expected "
                          << expected_flow << "), pressure "
                          << pres->GetTuple1(i) << " (expected "
                          << expected_pressure << ")" << std::endl;
            }
            num_failed++;
        }
    }
    if (num_failed)
    {
        std::cerr << "ERROR: " << num_failed << " of "
                  << output->GetNumberOfPoints()
                  << " points differ from the analytic solution" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "All " << output->GetNumberOfPoints()
              << " points match the analytic flow " << expected_flow
              << " and pressure" << std::endl;
    return EXIT_SUCCESS;
}

int compare(const std::string& reference_file, const std::string& output_file)
{
    auto reference = read_output(reference_file);
    auto output = read_output(output_file);
    if (!reference || !output) return EXIT_FAILURE;
    if (reference->GetNumberOfPoints() != output->GetNumberOfPoints())
    {
        std::cerr << "ERROR: number of points differs" << std::endl;
        return EXIT_FAILURE;
    }

    // Results must agree with the reference to rounding.
    int num_failed = 0;
    double max_rel_diff = 0.0;
    for (auto name : {"Velocity", "Pressure"})
    {
        auto ref = reference->GetPointData()->GetArray(name);
        auto out = output->GetPointData()->GetArray(name);
        for (vtkIdType i = 0; i < reference->GetNumberOfPoints(); i++)
        {
            double a[3], b[3];
            reference->GetPoint(i, a);
            output->GetPoint(i, b);
            bool same_point = a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
            double r = ref->GetTuple1(i);
            double o = out->GetTuple1(i);
            double rel_diff = std::abs(o - r) / std::max(std::abs(r), 1e-300);
            max_rel_diff = std::max(max_rel_diff, rel_diff);
            if (!same_point || !(rel_diff <= compare_tolerance))
            {
                if (num_failed < 10)
                {
                    std::cerr << std::setprecision(17) << name << " at point "
                              << i << ": " << r << " (reference) vs " << o
                              << " (relative difference " << rel_diff << ")"
                              << std::endl;
                }
                num_failed++;
            }
        }
    }
    if (num_failed)
    {
        std::cerr << "ERROR: '" << output_file << "' differs from '"
                  << reference_file << "' in " << num_failed
                  << " values (max relative difference " << max_rel_diff << ")"
                  << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "'" << output_file << "' matches '" << reference_file
              << "' (max relative difference " << max_rel_diff << ")"
              << std::endl;
    return EXIT_SUCCESS;
}
}  // namespace

int main(int argc, char* argv[])
{
    std::string command = argc > 1 ? argv[1] : "";
    if (command == "make-case" && argc == 4) return make_case(argv[2], argv[3]);
    if (command == "check-analytic" && argc == 3)
        return check_analytic(argv[2]);
    if (command == "compare" && argc == 4) return compare(argv[2], argv[3]);
    std::cerr << "Usage:\n"
              << "  svslicer_test_utils make-case RESULT.vtu CENTERLINE.vtp\n"
              << "  svslicer_test_utils check-analytic OUTPUT.vtp\n"
              << "  svslicer_test_utils compare REFERENCE.vtp OUTPUT.vtp"
              << std::endl;
    return EXIT_FAILURE;
}
