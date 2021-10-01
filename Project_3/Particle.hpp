#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <armadillo>

class Particle {
public:
    // member variables
    double q;               // charge of particle
    double m;               // mass of particle
    // constructor
    Particle(double q, double m);

    // member variables
    arma::vec r();            // position of particle
    arma::vec v();            // velocity of particle

};

#endif