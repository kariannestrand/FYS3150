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

    // making list/contatiner for particle objects

    for (int i = 0; i < n_; i++){
        particles_.push_back(Particle(q_vec(i), m_vec(i), R.col(i), V.col(i)));
    }

}

vec PenningTrap::external_B_field(){
    vec B = vec(3).fill(0.);
    B(2) = B0_;
    return B;
}


vec PenningTrap::external_E_field(int i){
    vec F = vec(3).fill(0);
    // position of E-field
    F(0) = -1.;
    F(1) = -1.;
    F(2) = 2.;

    vec r = particles_[i].r_;

    vec E = - V0_/(d_*d_)*F % r;      // E-field

    return E;
}

vec PenningTrap::force_particle(int i, int j){
    vec F;
    vec r = particles_[i].r_ - particles_[j].r_;

    double q_i = particles_[i].q_;
    double q_j = particles_[j].q_;
    vec dr = abs(r) % abs(r) % abs(r);

    F = ke_*(q_i*q_j)/dr % r;

    return F;
}

vec PenningTrap::total_force_external(int i){
    vec F = vec(3).fill(0);
    vec E = external_E_field(i);
    vec B = external_B_field();
    vec v = particles_[i].v_;

    double q = particles_[i].q_;

    F = q*E + cross(q*v, B);

    return F;
}

vec PenningTrap::total_force_particles(int i){
    vec F = vec(3).fill(0);
    for (int j = 0; j < n_; j++){
        if (i == j){
           continue;
        }
        F += force_particle(i, j);
    }

    return F;
}


vec PenningTrap::total_force(int i){
    vec F;
    F = total_force_particles(i) + total_force_external(i);

    return F;

}

/*
void PenningTrap::evolve_RK4(double dt){

}
*/

void PenningTrap::evolve_forward_Euler(double dt){

    // Euler

    for (int i = 0; i < n_; i++){
        Particle& p_i = particles_[i];
        F = total_force(i);
        
        p_i.r = p_i.r + p_i.v*dt;
        p_i.v = p_i.v + F/p_i.m*dt; 

        }
    }

}


