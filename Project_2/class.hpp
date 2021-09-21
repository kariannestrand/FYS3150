#ifndef CLASS_HPP
#define CLASS_HPP

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <assert.h>

// declear class or something

class MyClass {
private:
    /* declaration of variables only accessible from within the class*/
    int N_;
    double a_, d_; // Member variables, the _ has something to do with style and to prevent something
public:
    /* Declaration of constructor, destructor and class methods. */
    MyClass(int N, double a, double d); //Constructor
    //virtual ~MyClass(); //Destructor

    arma::mat num(); //Compute value y at point (x, y) given x.
};

#endif
