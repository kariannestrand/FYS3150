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

    int n = 2;                                  // number of particles
    int dim = 3;                                // dimension (x,y,z)

    double t = 100.;                            // total time, [mu*s]
    double dt = 0.001;                          // time step, [mu*s]
    int N = t/dt;                               // number of time steps

    vec q_vec = vec(n).fill(q);                 // vector with charges
    vec m_vec = vec(n).fill(m);                 // vector with masses

    mat pos = mat(dim, n).randn() - 0.5*d;      // fill in initial conditions for position here, just have random values for now
    mat vel = mat(dim, n).randn() - 0.5*d;      // fill in initial conditions for position here, just have random values for now

    double omega_0 = q*B0/m;
    double omega_z = sqrt(2*q*V0/(m*d*d));

    double v_0 = 1.;
    double x_0 = 1.;
    double omega_p = (1/2.)*(omega_0 + sqrt(omega_0*omega_0 - 2*omega_z*omega_z));
    double omega_m = (1/2.)*(omega_0 - sqrt(omega_0*omega_0 - 2*omega_z*omega_z));

    double Ap = (v_0 + omega_p*x_0)/(omega_p - omega_m);
    double Am = - (v_0 + omega_m*x_0)/(omega_p - omega_m);

    bool write = true;                          // creates txt-files if true
    bool interaction = false;                   // accounts for particle interactions if true
    bool euler = false;                         // runs evolve_forward_Euler method if true
    bool rk4 = true;                            // runs evolve_RK4 method if true
    bool analytical = false;                    // runs analytical_solution method if true


    if (euler){
        PenningTrap penningtrap0 = PenningTrap(B0, V0, d, ke, n, N, pos, vel, q_vec, m_vec, write, interaction);
        penningtrap0.evolve_forward_Euler(dt);
    }
    if (rk4){
        PenningTrap penningtrap1 = PenningTrap(B0, V0, d, ke, n, N, pos, vel, q_vec, m_vec, write, interaction);
        penningtrap1.evolve_RK4(dt);
    }
    if (analytical){
        PenningTrap penningtrap2 = PenningTrap(B0, V0, d, ke, n, N, pos, vel, q_vec, m_vec, write, interaction);
        penningtrap2.analytical_solution(dt, Ap, Am, omega_p, omega_m);
    }

    return 0;
}
