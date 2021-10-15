#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace arma;
using namespace std;

int main(int argc, char const *argv[]){
    double q = 1.;                              // charge of Ca+ particle, [e]
    double m = 40.078;                          // atomic mass of Ca+, [u]

    double B0 = 9.65e1;                         // magnetic field strength, [u/mu*s*e]
    double V0 = 9.65e8;                         // applied potential, [u*(mu*s)^2/(mu*s)^2*e]
    double d = 10e4;                            // characteristic dimension, [mu*m]

    double ke = 1.38935333e5;                   // Couloumb constant, [u*(mu*m)^3/(mu*s*e)^2]

    int n = 2;                    // number of particles
    int dim = 3;                  // dimension (x,y,z)

    double t = 100.;                            // total time, [mu*s]
    double dt = 0.001;                          // time step, [mu*s]
    int N = t/dt;                               // number of time steps

    bool write = true;                          // creates txt-files if true

    vec q_vec = vec(n).fill(q);                 // vector with charges
    vec m_vec = vec(n).fill(m);                 // vector with masses

    mat pos = mat(dim, n).randn() - 0.5*d;      // fill in initial conditions for position here, just have random values for now
    mat vel = mat(dim, n).randn() - 0.5*d;      // fill in initial conditions for position here, just have random values for now


    bool euler = false;
    if (euler){
        PenningTrap penningtrap0 = PenningTrap(B0, V0, d, ke, n, N, pos, vel, q_vec, m_vec);     // calling penningtrap
        penningtrap0.evolve_forward_Euler(dt, write);
    }

    bool rk4 = true;
    if (rk4){
        PenningTrap penningtrap1 = PenningTrap(B0, V0, d, ke, n, N, pos, vel, q_vec, m_vec);     // calling penningtrap
        penningtrap1.evolve_RK4(dt, write);
    }

    return 0;
}
