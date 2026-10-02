#include <iostream>
#include <fstream>
#include <string>
#include "src/solver.h"
#include "src/engine.h"

bool tempCases(double temp);

int main() {

    std::vector<double> temps;
    temps.clear();
    for (int i =0; i<11; i++){
        temps.push_back(1000+100*i);
    }
    for (int c =0; c<temps.size(); c++){
        tempCases(temps[c]);
    }


    return 0;
}

bool tempCases(double temp){

    // Constructing Inputs
    solver_inputs inputs;
    inputs.temperature = temp; //K
    inputs.surface_size_X = 220;
    inputs.surface_size_Y = 220;
    inputs.particle_flux = 50;
    inputs.time_step_size = 1e-6;
    inputs.time_step_no = 8500;

    // Initializing solver with inputs
    Solver solver(inputs);
    // Executing simulation with given inputs
    solver.execute();
 // Writing values to the txt files so we can post process data in MATLAB
    std::ofstream surfCovFile;
    std::string filename1 = "surfCovTot" + std::to_string((int) temp) +".txt";
    surfCovFile.open(filename1);
    for (double ti : solver.surface_cov){
        surfCovFile<<ti<<"\n";
    }
    surfCovFile.close();

    std::ofstream surfCovOFile;
    std::string filename2 = "surfCovO" + std::to_string((int) temp) +".txt";
    surfCovOFile.open(filename2);
    for (double ti : solver.surf_O){
        surfCovOFile<<ti<<"\n";
    }
    surfCovOFile.close();

    std::ofstream surfCovCOFile;
    std::string filename3 = "surfCovCO" + std::to_string((int) temp) +".txt";
    surfCovCOFile.open(filename3);
    for (double ti : solver.surf_CO){
        surfCovCOFile<<ti<<"\n";
    }
    surfCovCOFile.close();

    std::ofstream carbRemFile;
    std::string filename4 = "carbonFlux" + std::to_string((int) temp) +".txt";
    carbRemFile.open(filename4);
    for (double ti : solver.carbonflux){
        carbRemFile<<ti<<"\n";
    }
    carbRemFile.close();

    return true;

}
