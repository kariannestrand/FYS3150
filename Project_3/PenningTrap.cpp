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

    // making contatiner for particle objects

    for (int i = 0; i < n_; i++){
        particles_.push_back(Particle(q_vec(i), m_vec(i), R.col(i), V.col(i)));
    }

}

vec PenningTrap::external_B_field(vec r){
    vec B = vec(3).fill(0.);
    B(2) = B0_;
    return B;
}

<<<<<<< HEAD
vec PenningTrap::external_E_field(vec r){
    vec F = vec(3).fill(0.);
    // position of E-field
    F(0) = -1.;
=======
vec PenningTrap::external_E_field(int i){
    vec F = vec(3).fill(0.); 
    // position of E-field       
    F(0) = -1.;                     
>>>>>>> 0d928b2aeac6f7a1ce5159934a9ba928810924da
    F(1) = -1.;
    F(2) = 2.;
    
    vec E = (-V0_/(d_*d_)*F) % particles_[i].r_;      // E-field

    return E;
}

vec PenningTrap::force_particle(int i, int j){
    vec F = vec(3).fill(0.);
<<<<<<< HEAD
    for (int i = 0; i < n_; i++){
        for (int j = 0; j < n_; j++){
            if (i !=j){
            F = ke_ * (particles_[i].q_ * particles_[j].q_) / ((particles_[i].r_ % particles_[i].r_) + (particles_[j].r_ % particles_[j].r_)) % abs((particles_[i].r_ - particles_[j].r_));
            }
        }
    }

=======
    vec r = particles_[i].r_ - particles_[j].r_;
    F = ke_ * (particles_[i].q_ * particles_[j].q_) / (abs(r)%abs(r)%abs(r)) % r; 
>>>>>>> 0d928b2aeac6f7a1ce5159934a9ba928810924da
    return F;
}

vec PenningTrap::total_force_external(vec E, vec B, int i){
    vec F = vec(3).fill(0.);
<<<<<<< HEAD
    for (int i = 0; i < n_; i++){
        F = particles_[i].q_ * E + cross(particles_[i].q_ * particles_[i].v_, B);
    }
=======
    F = particles_[i].q_ * E + cross(particles_[i].q_ * particles_[i].v_, B);

>>>>>>> 0d928b2aeac6f7a1ce5159934a9ba928810924da
    return F;
}



<<<<<<< HEAD
}
=======
vec PenningTrap::total_force_particles(vec F_particles){
    vec F;
    for (int i = 0; i < n_; i++){
        for (int j = 0; j < n_; j++){
            if (i != j){
                F += F_particles;
            }
        }
    }

    return F;
}  
>>>>>>> 0d928b2aeac6f7a1ce5159934a9ba928810924da

/*
vec PenningTrap::total_force(int i){

}


void PenningTrap::evolve_RK4(double dt){

}

void PenningTrap::evolve_forward_Euler(double dt){

}
*/
