#ifndef FUNCTIONS_
#define FUNCTIONS_

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <complex>

//declaring functions
arma::cx_vec state(arma::cx_mat U_in);
arma::cx_mat matrix(double r, arma::cx_vec a, arma::cx_mat A);

#endif
