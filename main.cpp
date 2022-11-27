#include <chrono>
#include <iostream>
#include <string>

#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>

#include "mesh.hpp"

int main(int argc, char* argv[]) {

  auto start = std::chrono::high_resolution_clock::now();
  if (argc < 3) {
    std::cout << "Not enough input arguments." << std::endl;
    return EXIT_FAILURE;
  }

  std::string vol_file_name = argv[1];
  std::cout << "Reading volume file: " << vol_file_name << std::endl;
  auto mesh = new Mesh();
  mesh->read_mesh(vol_file_name);

  // Read in the centerline file.
  std::string cl_file_name = argv[2];
  std::cout << "Reading centerline file: " << cl_file_name << std::endl;
  auto cl_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  cl_reader->SetFileName(cl_file_name.c_str());
  cl_reader->Update();
  auto cl_polydata = vtkSmartPointer<vtkPolyData>::New();
  cl_polydata->DeepCopy(cl_reader->GetOutput());

  int num_cl_point_arrays =
      cl_polydata->GetPointData()->GetNumberOfArrays();
  std::vector<std::string> remove_data_names;
  for (int i = 0; i < num_cl_point_arrays; i++) {
    auto name =
        std::string(cl_polydata->GetPointData()->GetArrayName(i));
    if ((startsWith(name, "velocity") == true) || (startsWith(name, "pressure") == true)) {
      remove_data_names.push_back(name);
    }
  }
  for (auto const& name : remove_data_names) {
    std::cout << "Removing array from centerline: " << name << std::endl;
    cl_polydata->GetPointData()->RemoveArray(name.c_str());
  }

  std::cout << "Extracting slices: " << std::endl;
  mesh->map_on_centerline(cl_polydata.GetPointer(), true, false);

  std::string output_file = argv[3];
  std::cout << "Writing result to: " << output_file << std::endl;
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(output_file.c_str());
  writer->SetInputData(cl_polydata);
  writer->Update();
  writer->Write();


  delete mesh;
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = stop - start;
  std::cout << "Completed in " << duration.count()<< "s" << std::endl;
}
