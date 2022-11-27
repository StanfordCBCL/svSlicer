/**
 * @file main.cpp
 * @brief Main slicer routine.
 * 
 */
#include <iostream>
#include <filesystem>
#include <string>
#include <stdexcept>

#include <vtkXMLPolyDataReader.h>
#include <vtkCleanPolyData.h>

/**
 * @brief Check if given string starts with a certain prefix
 * 
 * @param string The string to check for the given prefix
 * @param prefix The prefix
 * @return true If string starts with the prefix
 * @return false If string does not start with the prefix
 */
bool starts_with(std::string& string, std::string prefix)
{
    if(string.find(prefix) == 0)
        return true;
    else
        return false;
}

/**
 * @brief 
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */
int main(int argc, char* argv[]) {

    if (argc < 3) {
        std::cerr << "ERROR: Not enough input arguments." << std::endl;
        return EXIT_FAILURE;
    }

    // Check if result file exists
    std::string result_file_name = argv[1];
    if (std::filesystem::exists(result_file_name) == false)
    {
        std::cerr << "ERROR: 3D result file '" << result_file_name << "' does not exist" << std::endl;
        return EXIT_FAILURE;
    }

    // Check if centerline file exists
    std::string cl_file_name = argv[2];
    if (std::filesystem::exists(cl_file_name) == false)
    {
        std::cerr << "ERROR: Centerline file '" << cl_file_name << "' does not exist" << std::endl;
        return EXIT_FAILURE;
    }

    // Read in the centerline file.
    std::cout << "Reading centerline file '" << cl_file_name << "'" << std::endl;
    auto cl_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    cl_reader->SetFileName(cl_file_name.c_str());
    cl_reader->Update();
    auto cl_polydata = vtkSmartPointer<vtkPolyData>::New();
    cl_polydata->DeepCopy(cl_reader->GetOutput());
    std::cout << cl_polydata->GetNumberOfVerts() << std::endl;
}
