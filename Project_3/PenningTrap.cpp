#include "PenningTrap.hpp"
#include "Particle.hpp"

using namespace arma;
using namespace std;

PenningTrap::PenningTrap(double B0, double V0, double d, double ke, int n, mat R, mat V, vec q_vec, vec m_vec){
    B0_ = B0;                   // magnetic field strength
    V0_ = V0;                   // applied potential
    d_ = d;                     // characteristic dimension
    ke_ = ke;
    n_ = n;

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
    vec F = vec(3).fill(0.);
    // position of E-field
    F(0) = -1.;
    F(1) = -1.;
    F(2) = 2.;
    vec E = (-V0_/(d_*d_)*F);      // E-field

    return E;
}

vec PenningTrap::force_particle(){
    vec F = vec(3).fill(0.);
    for (int i = 0; i < n_; i++){
        for (int j = 0; j < n_; j++){
            if (i !=j){
            F = ke_ * (particles_[i].q_ * particles_[j].q_) / ((particles_[i].r_ % particles_[i].r_) + (particles_[j].r_ % particles_[j].r_)) % abs((particles_[i].r_ - particles_[j].r_));
            }
        }
    }

    return F;
}

vec PenningTrap::total_force_external(){
    vec F = vec(3).fill(0.);
    for (int i = 0; i < n_; i++){
        F = particles_[i].q_ * E + cross(particles_[i].q_ * particles_[i].v_, B);
    }
    return F;
}


/*
vec PenningTrap::total_force_particles(int i){

}

vec PenningTrap::total_force(int i){

}

void PenningTrap::evolve_RK4(double dt){

}

void PenningTrap::evolve_forward_Euler(double dt){

}
*/
