//
// Created by Kaan Kirmanoglu on 11/22/21.
//

#ifndef PROJECTPROGRAM_SOLVER_H
#define PROJECTPROGRAM_SOLVER_H
#include <vector>
#include "solver_inputs.h"
#include "Constants.h"
#include <iostream>
#include <cmath>
#include "engine.h"
#include "particle.h"
#include <list>


class Solver {
public:
    Solver(solver_inputs inputs); // Initializing constructer
    std::vector<std::vector<int> > surface_sites; // Surface sites grid
    bool execute(); // Function that runs simulatiin
    // Data to be collected and processed
    std::vector<double> surface_cov;
    std::vector<double> surf_O;
    std::vector<double> surf_CO;
    std::vector<int> carbonflux;

private:

    bool applyGasCollisons(); // Perform gas collisions and gas surface reactions
    int gsReactid(); // Select which gas surface reactions to performed
    bool surfaceStep(); // Goes through adsorbed particles and desorbs when time counter exceeds desorption time
    bool deleteParticles();
    bool recordData();
    double temp; //temperature (K)
    double press;
    // Adsorbed number of total, O and CO in surface
    int adsCount;
    int adsCountO;
    int adsCountCO;
    double dt; // time step size
    double total_time; // number of time steps run
    int coll_per_step; //collisions per time step
    double stick_coeff = 1;
    double k_des=2;
    double k_desO;
    double k_desCO;
    std::vector<double> k_AMGS; //gas surface Arrhenius reaction rates [0]= LH3 O{a}, [1]= LH3 CO{a},
    // [2] = LH1 O and [3] = LH1 CO
    Engine rng; //random number generator
    int X;
    int Y;
    std::list<Particle> particles; // list of LH3 formed particles adsorbed, added and deleted as sim goes on
    double mSinv; //how many total site there is, adsCount*mSinverse gives surface coverage
    int cRemoved; //number of carbon atoms removed in a time step


};


#endif //PROJECTPROGRAM_SOLVER_H
