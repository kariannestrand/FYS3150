#ifndef FUNCTIONS_
#define FUNCTIONS_

#include <armadillo>
#include <iostream>
#include <iomanip>

// Declearing functions
arma::vec exact(arma::vec x);
void forward(arma::vec a, arma::vec *b, arma::vec c, arma::vec *g, int n);
void backward(arma::vec b, arma::vec c, arma::vec g, arma::vec *v, int n);
void gauss_elim_gen(arma::vec a, arma::vec b, arma::vec c, arma::vec g, arma::vec *v, int n);
void writetofile_exact(arma::vec x, arma::vec u, int n);
void writetofile_approx(arma::vec x, arma::vec v, int n);
void writetofile(arma::vec x, arma::vec u, arma::vec v, int n);

#endif
