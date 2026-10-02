//
// Created by Kaan Kirmanoglu on 11/22/21.
//

#ifndef PROJECTPROGRAM_PARTICLE_H
#define PROJECTPROGRAM_PARTICLE_H

// Surface Particles
class Particle {
public:
    int xLoc; // locations on site
    int yLoc;
    int type; // type 0 = O{a}, 1 = CO{a}
    double tau; // time counter
    double t_des; // desorption time
    bool skip; // if skip == true, particle desorped and it will be removed from list
};


#endif //PROJECTPROGRAM_PARTICLE_H
