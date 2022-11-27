#include <iostream>
#include <string>

#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>

int main(int argc, char* argv[]) {

    if (argc < 3) {
        std::cout << "Not enough input arguments." << std::endl;
        return EXIT_FAILURE;
    }

}