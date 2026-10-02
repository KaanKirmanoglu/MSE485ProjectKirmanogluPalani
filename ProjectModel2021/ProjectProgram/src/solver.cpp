//
// Created by Kaan Kirmanoglu on 11/22/21.
//

#include "solver.h"
Solver::Solver(solver_inputs inputs) {
        X = inputs.surface_size_X;
        Y = inputs.surface_size_Y;
        surface_sites.clear();
        surface_sites.resize(inputs.surface_size_X);
        double maxStick = X*Y;
        mSinv = 1/(maxStick);
        for (auto & surface_site : surface_sites){
            surface_site.resize(inputs.surface_size_Y);
        }
        for (int i=0;i<surface_sites.size(); i++){
            for (int j =0; j<surface_sites[i].size(); j++){
                surface_sites[i][j] =0;
            }
        }
        std::cout<<"Surface Temperature: "<<inputs.temperature<<std::endl;
        //std::cout<<"Pressure: "<< inputs.pressure<<std::endl;
        temp = inputs.temperature;
        press = inputs.pressure;
        coll_per_step = inputs.particle_flux;
        rng.seed(1);
        dt = inputs.time_step_size;
        total_time= inputs.time_step_no;
        surface_cov.clear();
        surface_cov.push_back(0);
        surf_CO.clear(); surf_CO.push_back(0);
        surf_O.clear(); surf_O.push_back(0);
        // Arrhenius rates taken from Swaminathan-Gopalan, et al. (2018). Development and validation of a
    // finite-rate model for carbon oxidation by atomic oxygen.
        k_desO = 0.05*temp*temp*exp(-3177.2/temp);
        k_desCO = 4485*exp(-1581.4/temp);
        k_AMGS.resize(5);
        k_AMGS[0]= 1; k_AMGS[1]= 153.0*exp(-4172.8/temp); k_AMGS[2] = 20.9*exp(-2449.3/temp);
        k_AMGS[3] = 1574.9*exp(-6240.0/temp); //k_AMGS[4] = 536.3*exp(-655.6/temp);
        adsCountO =0; adsCountCO =0; adsCount = 0; cRemoved = 0;
        carbonflux.clear(); carbonflux.push_back(0);

}

bool Solver::execute() {
    for (int ti =0; ti<total_time; ti++){
        surfaceStep(); // Surface step put first because otherwise time counter would add dt to particles
        // right after they adsorb
        applyGasCollisons();
        deleteParticles();
        recordData();
        std::cout<<"Timestep: "<<ti<<std::endl;
    }

    return true;
}
//Selecting random number on surface to represent gas collisions,
// checked if site is available then a reaction selected stochastically performed
bool Solver::applyGasCollisons() {
    //std::cout<<"Applying Gas collisions"<<std::endl;
    for (int p =0; p<coll_per_step; p++){
        int randX = floor((X-1)*rng.unit_uniform()+0.5);
        int randY = floor((Y-1)*rng.unit_uniform()+0.5);
       // std::cout<<"rndX = "<<randX<<" rndY = "<<randY<<std::endl;
        if (surface_sites[randX][randY] != 0){
            continue;
        }
        else{
            int reactid = gsReactid();
            if (reactid == 0 ){
                surface_sites[randX][randY] = 1;
                adsCount++;
                adsCountO++;
                Particle tempParticle;
                tempParticle.type = 0;
                tempParticle.xLoc = randX;
                tempParticle.yLoc = randY;
                tempParticle.tau = 0;
                tempParticle.t_des = -std::log(rng.unit_uniform())/k_desO;
                tempParticle.skip =false;
                particles.push_back(tempParticle);
            }
            if (reactid == 1 ){
                surface_sites[randX][randY] = 2;
                adsCount++;
                adsCountCO++;
                Particle tempParticle;
                tempParticle.type = 1;
                tempParticle.xLoc = randX;
                tempParticle.yLoc = randY;
                tempParticle.tau = 0;
                tempParticle.t_des = -std::log(rng.unit_uniform())/k_desCO;
                tempParticle.skip = false;
                particles.push_back(tempParticle);
            }
            if (reactid == 3 ){
                cRemoved++;
            }

        }
    }
    return true;
}


// desorbed particles marked for deletion, surface site is freed, adsorbed number of particles and
// carbon removal rates are adjusted
bool Solver::surfaceStep() {
    for (auto it = particles.begin(); it!=particles.end(); ++it){
        it->tau += dt;
        if (it->tau > it->t_des){
            it->skip = true;
            adsCount--;
            if (it->type == 1){
                adsCountCO--;
                cRemoved++;
            }
            else{
                adsCountO--;
            }
            surface_sites[it->xLoc][it->yLoc] = 0;
        }
    }
    return true;
}

bool Solver::deleteParticles() {
    std::list<Particle>::const_iterator itr = particles.cbegin();
    while (itr != particles.cend())
    {
        // remove all particles that have been skipped
        if (itr->skip) {
            itr = particles.erase(itr);
        }
        else {
            ++itr;
        }
    }
    return true;
}


// Selecting which reaction to be performed stochastically
int Solver::gsReactid() {
    double rand = rng.unit_uniform();
    std::vector<double> P_AMGS; P_AMGS.resize(5);
    double sumProb = 0;
    P_AMGS[0] = k_AMGS[0]; P_AMGS[1] = mSinv*adsCount*k_AMGS[1]; P_AMGS[2] = k_AMGS[2];
    P_AMGS[3] = mSinv*adsCount*k_AMGS[3];

    sumProb = P_AMGS[0] + P_AMGS[1] + P_AMGS[2] + P_AMGS[3];
    int id = -1;
    double probs =0;
    for (int rx =0 ; rx<4 ; rx++){
        double normProbi = P_AMGS[rx]/sumProb;
        if (rand<normProbi+probs){
            id = rx;
            break;
        }
        probs += normProbi;
    }
   // std::cout<<"reaction id: "<<id<<std::endl;
    return id;
}

bool Solver::recordData() {
    std::cout<<"adsCountTot: "<<adsCount<<" surfCov: "<<mSinv*adsCount <<std::endl;
    std::cout<<"adsCountO: "<<adsCountO<<" surfCov: "<<mSinv*adsCountO<<std::endl;
    std::cout<<"adsCountCO: "<<adsCountCO<<" surfCov: "<<mSinv*adsCountCO<<std::endl;
    std::cout<<"carbon rem: "<<cRemoved<<std::endl;
    surf_O.push_back(mSinv*adsCountO);
    surf_CO.push_back(mSinv*adsCountCO);
    surface_cov.push_back(mSinv*adsCount);
    carbonflux.push_back(cRemoved);
    cRemoved =0;
    return true;
}