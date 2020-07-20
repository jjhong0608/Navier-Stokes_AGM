// Axial Line Generator Hear files
#include "Header/AXLGEN.hpp"
#include <cstdlib>

int main(int argc, char const *argv[]) {
    int nIter;
    double initialTime;
    double terminalTime;
    double dt;
    int timeStep;
    string Geometry_filename, AGL_output_file;
    string AxialFile, AGM_output_file;
    string AGL_output_fileTEMP;
    string AGM_output_fileTEMP;
    string property;
    ifstream InputFile(argv[1]);

    if (argc != 2) {
        printf("%s\n", "Usage: ./AGM <Input file>");
        return 1;
    }

    if (!InputFile.is_open()) {
        printf("%s%s\n", "No Input file: ", argv[1]);
        printf("%s\n", "Please Check Input file name");
        return 1;
    }

    InputFile >> property >> Geometry_filename;
    InputFile >> property >> AGL_output_fileTEMP;
    InputFile >> property >> AGM_output_fileTEMP;
    InputFile >> property >> initialTime;
    InputFile >> property >> terminalTime;
    InputFile >> property >> timeStep;
    InputFile >> property >> dt;
    InputFile >> property >> nIter;

    printf("%s\n", "Run the Axial Line Generator");
    AXLGEN(Geometry_filename, AGL_output_fileTEMP, AGM_output_fileTEMP, initialTime, terminalTime, timeStep, dt, nIter);

    return 0;
}
