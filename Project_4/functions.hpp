#ifndef FUNCTIONS_
#define FUNCTIONS_

#include <armadillo>
#include <iostream>

// Declearing functions:
arma::mat spin_matrix(int L);
void initialize(int L, arma::mat &S, double &E, double &M, int N);
int delta_E(arma::mat &S, int L, int i, int j);
void metropolis(arma::mat &S, int L, double T, double &E, double &M, int N_cycles, int N);
inline int PBC(int i, int limit, int add);

#endif
