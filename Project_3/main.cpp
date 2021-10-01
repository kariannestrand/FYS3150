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

    int n = 3;                    // number of particles
    int dim = 3;                  // dimension (x,y,z)

    vec q_vec = vec(n).fill(q);     // vector with charges
    vec m_vec = vec(n).fill(m);     // vector with masses

    mat R = mat(dim, n).randn();        // fill in initial conditions for position here, just have random values for now
    mat V = mat(dim, n).randn();        // fill in initial conditions for position here, just have random values for now

    PenningTrap penningtrap = PenningTrap(B0, V0, d, R, V, q_vec, m_vec);    // calling penningtrap

// printing individual masses, charges, positions and velocities (we dont need this but nice to see what the particles-container contains)
    double x;
    double y;
    for (int i = 0; i < n; i++){
        cout << "position:" << endl;
        penningtrap.particles_[i].r_.print();
        cout << "velocity:" << endl;
        penningtrap.particles_[i].v_.print();
        x = penningtrap.particles_[i].q_;
        cout << "charge:" << x << endl;
        y = penningtrap.particles_[i].m_;
        cout << "mass:" << y << endl;
    }

    return 0;

}
