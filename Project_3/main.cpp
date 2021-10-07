#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace arma;
using namespace std;

int main(int argc, char const *argv[])
{
    double q = 20;                // charge of Ca+ particle [e]
    double m = 40.078;            // atomic mass of Ca+ [u]

    double B0 = 9.65e1;           // magnetic field strength 
    double V0 = 9.65e8;           // applied potential 
    double d = 10e4;              // characteristic dimension

    double ke = 1.38935333e5;     // Couloumb constant

    int n = 2;                    // number of particles
    int dim = 3;                  // dimension (x,y,z)

    vec r = vec(3).fill(0);      // initial condition for position (filled with zeros for now)

    vec q_vec = vec(n).fill(q);     // vector with charges
    vec m_vec = vec(n).fill(m);     // vector with masses

    mat R = mat(dim, n).randn();        // fill in initial conditions for position here, just have random values for now
    mat V = mat(dim, n).randn();        // fill in initial conditions for position here, just have random values for now

    PenningTrap penningtrap = PenningTrap(B0, V0, d, ke, n, R, V, q_vec, m_vec);    // calling penningtrap

    vec B = penningtrap.external_B_field(r);        // external B-field
    vec E = penningtrap.external_E_field(r);        // external E-field
    vec F_particles = penningtrap.force_particle();
    vec F_total_ext = penningtrap.total_force_external();
    F_total_ext.print();


    return 0;

}
