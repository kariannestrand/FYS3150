#include "PenningTrap.hpp"
#include "Particle.hpp"

using namespace arma;
using namespace std;

PenningTrap::PenningTrap(double B0, double V0, double d, mat R, mat V, vec q_vec, vec m_vec){
    B0_ = B0;                   // magnetic field strength 
    V0_ = V0;                   // applied potential 
    d_ = d;                     // characteristic dimension

    int n_particles = R.n_cols;         // number of particles is number of columns in matrix with positions

    // making contatiner for particle objects

    for (int i = 0; i < n_particles; i++){
        particles_.push_back(Particle(q_vec(i), m_vec(i), R.col(i), V.col(i)));
    }

}

vec PenningTrap::external_B_field(vec r){
    vec B = vec(3).fill(0.);
    B(2) = B0_;
    return B;
}

vec PenningTrap::external_E_field(vec r){
    vec F = vec(3).fill(0.);        // random vector to get E-field
    F(0) = -1.;                     // got these values from gradient of V
    F(1) = -1.;
    F(2) = 2.;
    vec E = (-V0_/(d_*d_)*F)%particles_[0].r_;      // element-wise multiplication with position-vector
    return E;
}