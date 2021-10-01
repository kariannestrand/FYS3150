#include "PenningTrap.hpp"

using namespace arma;
using namespace std;

PenningTrap::PenningTrap(double B0, double V0, double d, arma::mat R, arma::mat V, arma::vec q_vec, arma::vec m_vec){
    B0_ = B0;                   // magnetic field strength 
    V0_ = V0;                   // applied potential 
    d_ = d;                     // characteristic dimension

    int n_particles = R.n_cols;         // number of particles is number of columns in matrix with positions

    // making contatiner for particle objects

    for (int i = 0; i < n_particles; i++){
        particles_.push_back(Particle(q_vec(i), m_vec(i), R.col(i), V.col(i)));
    }

}