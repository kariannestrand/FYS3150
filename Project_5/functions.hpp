#ifndef FUNCTIONS_
#define FUNCTIONS_

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>

// declaring functions
arma::cx_mat initial(double mean_x, double mean_y, double var_x, double var_y, double p_x, double p_y, double M, double h);
arma::cx_mat potential(double v0, double M, int n_slits, double slit_size, double separation_size, double wall_thickness);
void vector_ab(arma::cx_double r, double dt, double M, arma::cx_vec &a, arma::cx_vec &b, arma::cx_mat V);
void matrix(arma::cx_double r, arma::cx_vec a, arma::cx_vec b, arma::sp_cx_mat &A, arma::sp_cx_mat &B, double M);
arma::cx_cube CrankNicolson(arma::cx_mat U, const arma::sp_cx_mat &B, const arma::sp_cx_mat &A, double M, int N);
void write_to_file(arma::cx_cube U_cube, std::string name);

#endif
