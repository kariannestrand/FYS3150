#ifndef FUNCTIONS_
#define FUNCTIONS_

#include <armadillo>
#include <iostream>
#include <iomanip>

// Declearing functions:
inline int PBC(int i, int limit, int add);
arma::mat spin_matrix(int L);
void initialize(int L, arma::mat &S, double &E, double &M, int N, bool random);
int delta_E(arma::mat &S, int L, int i, int j);
void metropolis(arma::mat &S, int L, double T, double &E, double &M, int N_cycles, int N, bool write, int burnin);

#endif
