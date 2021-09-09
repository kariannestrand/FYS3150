#ifndef FUNCTIONS_
#define FUNCTIONS_

#include <armadillo>
#include <iostream>
#include <iomanip>

// Declaring functions
arma::vec exact(arma::vec x);
void forward(arma::vec a, arma::vec *b, arma::vec c, arma::vec *, int);
void backward(arma::vec, arma::vec, arma::vec, arma::vec*, int);
void elimination(arma::vec, arma::vec, arma::vec, arma::vec, arma::vec*, int);
void writetofile_exact(arma::vec, arma::vec);
void writetofile_approx(arma::vec, arma::vec, int);

#endif