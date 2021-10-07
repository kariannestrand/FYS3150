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

    // making lis/contatiner for particle objects

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

void PenningTrap::evolve_forward_Euler(double dt){
    vec a;
    vec v;
    vec r;

    double omega_p = ;
    double omega_m = ;

    double A_p = ;
    double A_m = ;

    double x0 = A_p + A_m;
    double v0 = ;

    v(0) = (0, v0, 0);
    r(0) = (x0, 0, z0);

    a(i) = F(i)/m;
    v(i+1) = v(i) + a(i)*dt;
    r(i+1) = r(i) + v(i+1)*dt;

    return r;

}

*/
