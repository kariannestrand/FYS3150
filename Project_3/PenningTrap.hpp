#ifndef PENNINGTRAP_HPP
#define PENNINGTRAP_HPP
  
#include
  
class PenningTrap {
public:
    // member variables
    double B0_in;
    double V0_in;
    double d_in;
    // constructor
    PenningTrap(double B0_in, double V0_in, double d_in);
    // class methods

    void add_particle(Particle p_in);           // add a particle to the trap
    arma::vec external_E_field(arma::vec r);    // external electric field at point r=(x,y,z)  
    arma::vec external_B_field(arma::vec r);    // external magnetic field at point r=(x,y,z) 
    arma::vec force_particle(int i, int j);     // force on particle_i from particle_j
    arma::vec total_force_external(int i);      // the total force on particle_i from the external fields
    arma::vec total_force_particles(int i);     // the total force on particle_i from the other particles
    arma::vec total_force(int i);               // the total force on particle_i from both external fields and other particles
    void evolve_RK4(double dt);                 // evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_forward_Euler(double dt);       // evolve the system one time step (dt) using Forward Euler
};

#endif