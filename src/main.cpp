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
        std::cout << "Not enough input arguments." << std::endl;
        return EXIT_FAILURE;
    }

    // Check if centerline files exist
    std::string cl_file_name = argv[2];
    if (std::filesystem::exists(cl_file_name) == false)
    {
        throw std::runtime_error("Centerline file does not exist!");
    }

    // Read in the centerline file.
    std::cout << "Reading centerline file: " << cl_file_name << std::endl;
    auto cl_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    cl_reader->SetFileName(cl_file_name.c_str());
    cl_reader->Update();
    auto cl_polydata = vtkSmartPointer<vtkPolyData>::New();
    cl_polydata->DeepCopy(cl_reader->GetOutput());
}
