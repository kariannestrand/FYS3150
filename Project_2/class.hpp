#ifndef CLASS_HPP
#define CLASS_HPP

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <assert.h>

class MyClass {
private:
    int N_;
    double a_, d_;
public:
    MyClass(int N, double a, double d);

    arma::mat tridiag_matrix();
    arma::mat eigen_vectors();
    arma::vec eigen_values();
    double max_offdiag_symmetric(arma::mat &A, int &k, int &l);
    void rotation(arma::mat &A, arma::mat &R, int k, int l);
};

#endif
