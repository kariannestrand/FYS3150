#include "ising.hpp"

using namespace arma;
using namespace std;

mat SpinConfiguration(int L);
double MeanEnergy(int L, mat S);
double Magnetization(int L, mat S);

int main(int argc, char const *argv[]){
    int L = 2;                          // lattice length
    double N = L*L;                     // number of spins


    mat S = SpinConfiguration(L);
    double E = MeanEnergy(L, S);
    double M = Magnetization(L, S);

    cout << S << endl;
    cout << E << endl;
    cout << M << endl;


    return 0;
}

mat SpinConfiguration(int L){
    arma_rng::set_seed_random();        // sets seed for randu
    mat S = mat(L, L, fill::randu);     // initialize the lattice spin values

    for(int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++){
            if (S(i, j) >= 0.5){
                S(i, j) = 1.0;          // spin up
            }
            else{
                S(i, j) = - 1.0;        // spin down
            }
        }
    }
    return S;
}

double MeanEnergy(int L, mat S){
    double E = 0.;                      // initialize mean energy

    for(int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++){
            if (i != j){
                E -= S(i, j);
            }
        }
    }
    return E;
}

double Magnetization(int L, mat S){
    double M = 0.;                      // initialize magnetization

    for(int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++){
            if (i != j){
                M += S(i, j);
            }
        }
    }
    return M;
}
