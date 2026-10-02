//
// Created by Kaan Kirmanoglu on 11/22/21.
//

#ifndef PROJECTPROGRAM_SOLVER_INPUTS_H
#define PROJECTPROGRAM_SOLVER_INPUTS_H

struct solver_inputs{
    double temperature;
    double pressure;
    double gas_mass;
    int surface_size_X;
    int surface_size_Y;
    double time_step_size;
    int time_step_no ;
    int particle_flux;
};


#endif //PROJECTPROGRAM_SOLVER_INPUTS_H
