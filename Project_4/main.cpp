#include "functions.hpp"

using namespace arma;
using namespace std;


int main(int argc, char const *argv[]){
    int L = 2;                         // lattice length
    double N = L*L;                     // number of spins

    //arma_rng::set_seed_random();        // sets seed for randu
    //mat S = mat(L, L, fill::randu);     // initialize the lattice spin values with random values from 0 to 1
    //mat S = mat(L, L);
    mat S = zeros<mat>(L,L);
    double E = 0.;                      // initialize energy
    double M = 0.;                      // initialize magnetization

    Initialize(L, &S, &E, &M);

    double epsilon = E/N;               // energy per spin
    double m = M/N;                     // magnetization per spin

    return 0;
}

