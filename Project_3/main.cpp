#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace arma;
using namespace std;

int main(int argc, char const *argv[])
{
    double q = 1;                // charge of Ca+ particle [e]
    double m = 40.078;            // atomic mass of Ca+ [u]

    double B0 = 9.65e1;           // magnetic field strength
    double V0 = 9.65e8;           // applied potential
    double d = 10e4;              // characteristic dimension

    double ke = 1.38935333e5;     // Couloumb constant

    int n = 2;                    // number of particles
    int dim = 3;                  // dimension (x,y,z)

    int t = 10;
    int N = 100;
    double dt = t*(1./N);


    vec r = vec(3).fill(0);      // initial condition for position (filled with zeros for now)

    vec q_vec = vec(n).fill(q);     // vector with charges
    vec m_vec = vec(n).fill(m);     // vector with masses

    mat R = mat(dim, n).randn() - 1/2*d;        // fill in initial conditions for position here, just have random values for now
    mat V = mat(dim, n).randn() - 1/2*d;        // fill in initial conditions for position here, just have random values for now

    PenningTrap penningtrap = PenningTrap(B0, V0, d, ke, n, N, R, V, q_vec, m_vec);    // calling penningtrap
    
    bool write = true;
    penningtrap.evolve_forward_Euler(dt, write);

    /*
    bool write = true;
    //for (int i = 0; i < N; i++){
    penningtrap.evolve_forward_Euler(dt, write);

    /*
    for (int i = 0; i < N; i++){
        penningtrap.evolve_RK4(dt, write);

    }
    */

    return 0;

}
