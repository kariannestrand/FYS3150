#ifndef FUNCTIONS_
#define FUNCTIONS_

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <complex>

//declaring functions
// arma::cx_vec state(arma::cx_mat U_in);
void vector_ab(double r, double dt, int M, arma::cx_vec &a, arma::cx_vec &b);
void matrix(double r, arma::cx_vec a, arma::cx_vec b, arma::sp_cx_mat &A, arma::cx_mat &B, int M);
void solver(arma::cx_mat U_in, arma::cx_mat B, arma::sp_cx_mat A);

#endif
