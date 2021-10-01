#include "Particle.hpp"

using namespace arma;
using namespace std;

Particle::Particle(double q, double m, vec r, vec v){
    q_ = q;                 // charge of Ca+ particle
    m_ = m;                 // atomic mass of Ca+ 
    r_.swap(r);             // position
    v_.swap(v);             // velocity
}


