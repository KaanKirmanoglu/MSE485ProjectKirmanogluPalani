//
// Created by Kaan Kirmanoglu on 11/22/21.
//

#ifndef PROJECTPROGRAM_SOLVER_INPUTS_H
#define PROJECTPROGRAM_SOLVER_INPUTS_H

// All inputs have defaults so that no field is ever read uninitialized.
struct solver_inputs{
    double temperature = 1200.0;    // surface temperature (K)
    double pressure = 0.0;          // not used by the current model
    double gas_mass = 0.0;          // not used by the current model
    int surface_size_X = 100;       // lattice sites in x
    int surface_size_Y = 100;       // lattice sites in y
    double time_step_size = 1e-6;   // time step (s)
    int time_step_no = 1000;        // number of time steps
    int particle_flux = 40;         // O atoms hitting the surface per time step
    unsigned int seed = 1;          // random-number seed
    bool verbose = false;           // print coverage at every time step

    // Verification mode: pure Langmuir adsorption-desorption.
    // Every collision on a free site adsorbs (sticking coefficient 1) and
    // adsorbates desorb with constant rate langmuir_k_des (1/s).
    bool langmuir_test = false;
    double langmuir_k_des = 2.0;
};


#endif //PROJECTPROGRAM_SOLVER_INPUTS_H
