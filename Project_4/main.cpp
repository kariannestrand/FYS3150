#include "ising.hpp"

using namespace arma;
using namespace std;

/* WITHOUT POINTERS */
/*
mat SpinConfiguration(int L);
double MeanEnergy(int L, mat S);
double Magnetization(int L, mat S);
*/


/* WITH POINTERS */
void SEM(int L, mat *S, double *E, double *M);

int main(int argc, char const *argv[]){
    int L = 2;                          // lattice length
    double N = L*L;                     // number of spins


    /* WITHOUT POINTERS */
    /*
    mat S = SpinConfiguration(L);       // declearing spin configuration
    double E = MeanEnergy(L, S);        // declearing mean energy
    double M = Magnetization(L, S);     // declearing magnetization
    */


    /* WITH POINTERS */
    // /*
    arma_rng::set_seed_random();        // sets seed for randu
    mat S = mat(L, L, fill::randu);     // initialize the lattice spin values with random values from 0 to 1
    double E = 0.;                      // initialize mean energy
    double M = 0.;                      // initialize magnetization

    SEM(L, &S, &E, &M);
    // */

    cout << S << endl;
    cout << E << endl;
    cout << M << endl;

    return 0;
}


/* WITHOUT POINTERS */
/*
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
*/



/* WITH POINTERS */
// /*
void SEM(int L, mat *S, double *E, double *M){
    for(int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++){

            // spin configuration
            if ((*S)(i, j) >= 0.5){
                (*S)(i, j) = 1.0;          // spin up
            }
            else{
                (*S)(i, j) = - 1.0;        // spin down
            }

            // mean energy and magnetization
            if (i != j){
                (*E) -= (*S)(i, j);
                (*M) += (*S)(i, j);
            }
        }
    }
}
// */
