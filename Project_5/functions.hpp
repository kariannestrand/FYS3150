#ifndef FUNCTIONS_
#define FUNCTIONS_

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>

//declaring functions
arma::cx_mat initial(double mean_x, double mean_y, double var_x, double var_y, int M);
void vector_ab(double r, double dt, int M, arma::cx_vec &a, arma::cx_vec &b);
void matrix(double r, arma::cx_vec a, arma::cx_vec b, arma::sp_cx_mat &A, arma::cx_mat &B, int M);
arma::cx_vec solver(arma::cx_mat U_in, arma::cx_mat B, arma::sp_cx_mat A);

#endif
