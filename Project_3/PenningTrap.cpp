#include "PenningTrap.hpp"

using namespace arma;
using namespace std;

PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in){
    B0_in = 9.65e1         // magnetic field strength [T]
    V0_in = 9.65e10        // applied potential
    d_in = 10e4            // characteristic dimension
}