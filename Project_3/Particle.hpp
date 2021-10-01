#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <armadillo>
#include <iostream>

class Particle {
public:
    // member variables
    double q_;               // charge of particle
    double m_;               // mass of particle
    arma::vec r_;            // position of particle
    arma::vec v_;            // velocity of particle

    // constructor
    Particle(double q, double m, arma::vec r, arma::vec v);

};

#endif