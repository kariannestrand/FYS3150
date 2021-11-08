#ifndef FUNCTIONS_
#define FUNCTIONS_

#include <armadillo>
#include <iostream>

// Declearing functions:

void Initialize(int L, arma::mat *S, double *E, double *M);
inline int PBC(int i, int limit, int add);

#endif
