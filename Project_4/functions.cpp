#include "functions.hpp"

using namespace arma;
using namespace std;


// function to initialize spin configuration, energy and magnetization
void Initialize(int L, mat *S, double *E, double *M){
    for(int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){

            if ((*S)(i, j) >= 0.5){
                (*S)(i, j) = 1.0;          // spin configuration, spin up
            }
            else{
                (*S)(i, j) = - 1.0;        // spin configuration, spin down
            }

            (*E) -= (*S)(i, j) * ((*S)(PBC(i, L, -1), j) + (*S)(i, PBC(j, L, -1)));
            (*M) += (*S)(i, j);

        }
    }
}


// inline function for Periodic Boundary Conditions
inline int PBC(int i, int limit, int add){
    // metropolis sampling is used, periodic boundary conditions
    // i: base index, limit: highest legal index (L), add: number added to i
    // a % b gives the remainder of a/b
    return (i + limit + add) % (limit);
}



/* COMMENT FROM KARIANNE */
// Other people from last year flip spins and stuff, which i guess is part of the Metropolis-thing
// HOWEVER, I dont know if this is a part of our project or not?
