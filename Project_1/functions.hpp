#ifndef FUNCTIONS_
#define FUNCTIONS_

#include <armadillo>
#include <iostream>
#include <iomanip>

// Declearing functions
arma::vec exact(arma::vec x);
void forward(arma::vec a, arma::vec *b, arma::vec c, arma::vec *g, int n);
void backward(arma::vec b, arma::vec c, arma::vec g, arma::vec *v, int n);
void elimination(arma::vec a, arma::vec b, arma::vec c, arma::vec g, arma::vec *v, int n);
arma::vec abs_err(arma::vec v, arma::vec u);
arma::vec rel_err(arma::vec v, arma::vec u);
void writetofile_exact(arma::vec x, arma::vec u, int n);
void writetofile_approx(arma::vec x, arma::vec v, int n);
void writetofile_abs_err(arma::vec x, arma::vec Delta, int n);
void writetofile_rel_err(arma::vec x, arma::vec epsilon, int n);
void writetofile(arma::vec x, arma::vec u, arma::vec v, arma::vec Delta, arma::vec epsilon, int n);

#endif
