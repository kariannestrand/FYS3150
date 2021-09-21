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

    arma::mat num();
    arma::mat eigen_vectors();
};

#endif
