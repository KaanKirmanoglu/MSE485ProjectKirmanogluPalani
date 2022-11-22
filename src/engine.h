//
// Created by Kaan Kirmanoglu on 11/30/21.
//

#include "prng_engine.h"
#ifndef PROJECTPROGRAM_ENGINE_H
#define PROJECTPROGRAM_ENGINE_H
class Engine{
private:
    sitmo::prng_engine engine;

public:
    void seed(int random_seed){
        engine.seed(random_seed);
    }
    double unit_uniform() {
        return std::max((double)engine()/sitmo::prng_engine::max(),1e-9);
    }
    double uniform_minusone_one(){
        return unit_uniform() * 2. - 1.;
    }
    double double_range(double range) {
        return unit_uniform() * range;
    }



};


#endif //PROJECTPROGRAM_ENGINE_H
