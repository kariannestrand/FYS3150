#ifndef CLASS_HPP
#define CLASS_HPP

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <chrono>

// declear class or something

class MyClass {
private:
    /* declaration of variables only accessible from within the class*/
    double c0_, c1_; // Member variables, the _ has something to do with style and to prevent something
public:
    /* Declaration of constructor, destructor and class methods. */
    MyClass(double c0, double c1); //Constructor
    virtual ~MyClass(); //Destructor

    double compute_val(double x); //Compute value y at point (x, y) given x.
};


#endif
