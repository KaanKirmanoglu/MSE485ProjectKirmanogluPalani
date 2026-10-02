#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include "src/solver.h"
#include "src/engine.h"

template <typename T>
static void writeSeries(const std::string& filename, const std::vector<T>& data) {
    std::ofstream file(filename);
    for (const T& value : data) {
        file << value << "\n";
    }
}

static bool tempCase(double temp, int particle_flux, int time_steps);

// Usage: ProjectProgram [particle_flux] [time_steps]
// Sweeps T = 1000-2000 K and writes one set of output files per temperature.
int main(int argc, char* argv[]) {
    const int particle_flux = (argc > 1) ? std::atoi(argv[1]) : 50;
    const int time_steps    = (argc > 2) ? std::atoi(argv[2]) : 8500;

    for (int i = 0; i < 11; i++){
        tempCase(1000.0 + 100.0*i, particle_flux, time_steps);
    }

    return 0;
}

static bool tempCase(double temp, int particle_flux, int time_steps){

    // Constructing inputs
    solver_inputs inputs;
    inputs.temperature = temp; // K
    inputs.surface_size_X = 220;
    inputs.surface_size_Y = 220;
    inputs.particle_flux = particle_flux; // O atoms per time step
    inputs.time_step_size = 1e-6;         // s
    inputs.time_step_no = time_steps;

    // Initializing solver with inputs and running the simulation
    Solver solver(inputs);
    solver.execute();

    // One value per time step, read by postprocess.py
    const std::string T = std::to_string(static_cast<int>(temp));
    writeSeries("surfCovTot" + T + ".txt", solver.surface_cov);
    writeSeries("surfCovO"   + T + ".txt", solver.surf_O);
    writeSeries("surfCovCO"  + T + ".txt", solver.surf_CO);
    writeSeries("carbonFlux" + T + ".txt", solver.carbonflux);

    return true;
}
