#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace arma;
using namespace std;

int main(int argc, char const *argv[]){
    bool write = false;                           // creates txt-files if true
    bool interaction = false;                      // accounts for particle interactions if true
    bool modified = true;                         // runs the program with time-dependent electrical field if true (time-independent if false)
    bool euler = false;                           // runs evolve_forward_Euler method if true
    bool rk4 = true;                              // runs evolve_RK4 method if true

    int n = 100;                                  // number of particles

    double t = 500.;                              // total time, [mu*s]
    double dt = 0.008;                              // time step, [mu*s]
    int N = t/dt;                                 // number of time steps

    double V0 = (2.5e-4)*9.65e8;                  // applied potential, [u*(mu*m)^2/(mu*s)^2*e]
    double d = 0.05*1.0e4;                        // characteristic dimension, [mu*m]

    double f = 0.1;                               // amplitude
    double No = 115;                              // number of steps in omega_v vector
    vec omega_v = linspace(0.2, 2.5, No);         // angular frequency, [MHz]

    double q = 1.;                                // charge of Ca+ particle, [e]
    double m = 40.078;                            // atomic mass of Ca+, [u]
    double B0 = 9.65e1;                           // magnetic field strength, [u/mu*s*e]
    double ke = 1.38935333e5;                     // Couloumb constant, [u*(mu*m)^3/(mu*s*e)^2]
    int dim = 3;                                  // dimension (x,y,z)

    vec q_vec = vec(n).fill(q);                   // vector with charges
    vec m_vec = vec(n).fill(m);                   // vector with masses

    arma_rng::set_seed_random();                  // filling initial conditions with random values

    mat pos = mat(dim, n).randn()*0.1*d;          // initial conditions for position
    mat vel = mat(dim, n).randn()*0.1*d;          // initial conditions for velocity


    if (modified){
        for (int k = 0; k < omega_v.size(); k++){
            if (euler){
                PenningTrap penningtrap_euler = PenningTrap(B0, V0, d, ke, f, omega_v, n, N, pos, vel, q_vec, m_vec, write, interaction, modified);
                penningtrap_euler.evolve_forward_Euler(dt, k);
            }
            if (rk4){
                PenningTrap penningtrap_rk4 = PenningTrap(B0, V0, d, ke, f, omega_v, n, N, pos, vel, q_vec, m_vec, write, interaction, modified);
                penningtrap_rk4.evolve_RK4(dt, k);
            }
        }
    }
    else{
        int k = 0;
        if (euler){
            PenningTrap penningtrap_euler = PenningTrap(B0, V0, d, ke, f, omega_v, n, N, pos, vel, q_vec, m_vec, write, interaction, modified);
            penningtrap_euler.evolve_forward_Euler(dt, k);
        }
        if (rk4){
            PenningTrap penningtrap_rk4 = PenningTrap(B0, V0, d, ke, f, omega_v, n, N, pos, vel, q_vec, m_vec, write, interaction, modified);
            penningtrap_rk4.evolve_RK4(dt, k);
        }
    }

    return 0;
}
