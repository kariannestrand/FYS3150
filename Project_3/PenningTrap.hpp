#ifndef PENNINGTRAP_HPP
#define PENNINGTRAP_HPP

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include "Particle.hpp"

class PenningTrap {
public:
    // member variables
    double B0_;           // magnetic field strength
    double V0_;           // applied potential
    double d_;            // characteristic dimension
    double ke_;           // coloumbs constant
    double f_;            // amplitude
    arma::vec omega_v_;
    int n_;               // number of particles
    double N_;            // number of time steps
    bool write_;
    bool interaction_;
    bool modified_;
    std::vector<Particle> particles_;       // to contain all the Particle objects in the Penning trap

    // constructor
    PenningTrap(double B0, double V0, double d, double ke, double f, arma::vec omega_v, int n, double N, arma::vec r1, arma::vec v1, arma::vec r2, arma::vec v2, arma::mat pos, arma::mat vel, arma::vec q_vec, arma::vec m_vec, bool write, bool interaction, bool modified);
    
    // class methods
    arma::vec external_E_field(int i, int k, double dt);       // external electric field at point r=(x,y,z)
    arma::vec external_B_field(int i);                         // external magnetic field at point r=(x,y,z)
    arma::vec force_particle(int i, int j);                    // force on particle_i from particle_j
    arma::vec total_force_external(int i, int k, double t);    // the total force on particle_i from the external fields
    arma::vec total_force_particles(int i);                    // the total force on particle_i from the other particles
    arma::vec total_force(int i, int k, double dt);            // the total force on particle_i from both external fields and other particles
    int particles_trapped();                                   // number of particles in trap 
    void evolve_RK4(double dt, int k);                         // evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_forward_Euler(double dt, int k);               // evolve the system one time step (dt) using Forward Euler
    
};

#endif
