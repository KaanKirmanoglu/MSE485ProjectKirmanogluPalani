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
    explicit Solver(const solver_inputs& inputs); // Initializing constructor
    std::vector<std::vector<int> > surface_sites; // Surface sites grid (0 = empty, 1 = O(a), 2 = CO(a))
    bool execute(); // Function that runs the simulation
    // Data collected at every time step
    std::vector<double> surface_cov;
    std::vector<double> surf_O;
    std::vector<double> surf_CO;
    std::vector<int> carbonflux;

private:

    bool applyGasCollisions(); // Perform gas collisions and gas-surface reactions
    int gsReactid(); // Select which gas-surface reaction is performed
    bool surfaceStep(); // Goes through adsorbed particles and desorbs when time counter exceeds desorption time
    bool deleteParticles();
    bool recordData();
    void adsorb(int x, int y, int type); // Place an O(a) (type 0) or CO(a) (type 1) on site (x, y)
    double temp; // temperature (K)
    // Adsorbed number of total, O and CO on the surface
    int adsCount;
    int adsCountO;
    int adsCountCO;
    double dt; // time step size
    int total_time; // number of time steps run
    int coll_per_step; // collisions per time step
    bool verbose;
    bool langmuir_test;
    double k_desO;
    double k_desCO;
    std::vector<double> k_AMGS; // gas-surface Arrhenius reaction rates [0] = LH3 O(a), [1] = LH3 CO(a),
    // [2] = LH1 O and [3] = LH1 CO
    Engine rng; // random number generator
    int X;
    int Y;
    std::list<Particle> particles; // list of LH3-formed adsorbed particles, added and deleted as the sim goes on
    double mSinv; // 1 / total number of sites; adsCount*mSinv gives surface coverage
    int cRemoved; // number of carbon atoms removed in a time step


};


#endif //PROJECTPROGRAM_SOLVER_H
