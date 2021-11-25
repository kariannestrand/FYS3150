#ifndef FUNCTIONS_
#define FUNCTIONS_

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <complex>

//declaring functions
arma::cx_vec state(arma::cx_mat U_in);
void matrix(double r, arma::cx_vec a, arma::cx_vec b, arma::cx_mat &A, arma::cx_mat &B, int M);

#endif
