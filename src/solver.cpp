//
// Created by Kaan Kirmanoglu on 11/22/21.
//

#include "solver.h"
#include <algorithm>

Solver::Solver(const solver_inputs& inputs) {
        X = inputs.surface_size_X;
        Y = inputs.surface_size_Y;
        surface_sites.assign(X, std::vector<int>(Y, 0));
        mSinv = 1.0/(static_cast<double>(X)*Y);
        std::cout<<"Surface Temperature: "<<inputs.temperature<<" K"<<std::endl;
        temp = inputs.temperature;
        coll_per_step = inputs.particle_flux;
        verbose = inputs.verbose;
        langmuir_test = inputs.langmuir_test;
        rng.seed(inputs.seed);
        dt = inputs.time_step_size;
        total_time = inputs.time_step_no;
        surface_cov.assign(1, 0.0);
        surf_CO.assign(1, 0.0);
        surf_O.assign(1, 0.0);
        carbonflux.assign(1, 0);
        // Arrhenius rates taken from Swaminathan-Gopalan, et al. (2018). Development and validation of a
        // finite-rate model for carbon oxidation by atomic oxygen.
        k_desO = 0.05*temp*temp*std::exp(-3177.2/temp);
        k_desCO = 4485.5*std::exp(-1581.4/temp);
        k_AMGS.resize(4);
        k_AMGS[0] = 1;
        k_AMGS[1] = 153.0*std::exp(-4172.8/temp);
        k_AMGS[2] = 20.9*std::exp(-2449.3/temp);
        k_AMGS[3] = 1574.9*std::exp(-6240.0/temp);
        if (langmuir_test) {
            // Verification case: only adsorption of O(a), constant desorption rate
            k_desO = inputs.langmuir_k_des;
        }
        adsCountO = 0; adsCountCO = 0; adsCount = 0; cRemoved = 0;

}

bool Solver::execute() {
    for (int ti = 0; ti < total_time; ti++){
        surfaceStep(); // Surface step put first because otherwise the time counter would add dt to particles
        // right after they adsorb
        applyGasCollisions();
        deleteParticles();
        recordData();
        if (verbose) {
            std::cout<<"Timestep: "<<ti<<"  surfCov: "<<surface_cov.back()
                     <<"  O: "<<surf_O.back()<<"  CO: "<<surf_CO.back()
                     <<"  carbon removed: "<<carbonflux.back()<<"\n";
        }
    }

    return true;
}

void Solver::adsorb(int x, int y, int type) {
    surface_sites[x][y] = type + 1;
    adsCount++;
    if (type == 0) { adsCountO++; } else { adsCountCO++; }
    Particle p;
    p.type = type;
    p.xLoc = x;
    p.yLoc = y;
    p.tau = 0;
    p.t_des = -std::log(rng.unit_uniform())/(type == 0 ? k_desO : k_desCO);
    p.skip = false;
    particles.push_back(p);
}

// Selecting a random site on the surface to represent each gas collision.
// If the site is available, a reaction is selected stochastically and performed.
bool Solver::applyGasCollisions() {
    for (int p = 0; p < coll_per_step; p++){
        // Uniform site selection: every site, including edge sites, has probability 1/(X*Y)
        int randX = std::min(static_cast<int>(X*rng.unit_uniform()), X - 1);
        int randY = std::min(static_cast<int>(Y*rng.unit_uniform()), Y - 1);
        if (surface_sites[randX][randY] != 0){
            continue;
        }
        int reactid = langmuir_test ? 0 : gsReactid();
        if (reactid == 0){        // LH3 O(a) formation
            adsorb(randX, randY, 0);
        }
        else if (reactid == 1){   // LH3 CO(a) formation
            adsorb(randX, randY, 1);
        }
        else if (reactid == 3){   // LH1 CO formation and immediate desorption
            cRemoved++;
        }
        // reactid == 2: LH1 O formation/desorption, O leaves the surface, no carbon removed
    }
    return true;
}


// Desorbed particles are marked for deletion, their surface site is freed, and the adsorbed counts and
// carbon removal are updated
bool Solver::surfaceStep() {
    for (auto& particle : particles){
        particle.tau += dt;
        if (particle.tau > particle.t_des){
            particle.skip = true;
            adsCount--;
            if (particle.type == 1){
                adsCountCO--;
                cRemoved++;
            }
            else{
                adsCountO--;
            }
            surface_sites[particle.xLoc][particle.yLoc] = 0;
        }
    }
    return true;
}

bool Solver::deleteParticles() {
    // remove all particles that have desorbed
    particles.remove_if([](const Particle& p){ return p.skip; });
    return true;
}


// Selecting which reaction is performed stochastically, with probability proportional to its rate
int Solver::gsReactid() {
    const double theta = mSinv*adsCount;
    double P_AMGS[4];
    P_AMGS[0] = k_AMGS[0];
    P_AMGS[1] = theta*k_AMGS[1];
    P_AMGS[2] = k_AMGS[2];
    P_AMGS[3] = theta*k_AMGS[3];
    const double sumProb = P_AMGS[0] + P_AMGS[1] + P_AMGS[2] + P_AMGS[3];

    const double target = rng.unit_uniform()*sumProb;
    double cumulative = 0;
    for (int rx = 0; rx < 4; rx++){
        cumulative += P_AMGS[rx];
        if (target < cumulative){
            return rx;
        }
    }
    return 3; // guards against round-off when target == sumProb
}

bool Solver::recordData() {
    surf_O.push_back(mSinv*adsCountO);
    surf_CO.push_back(mSinv*adsCountCO);
    surface_cov.push_back(mSinv*adsCount);
    carbonflux.push_back(cRemoved);
    cRemoved = 0;
    return true;
}
