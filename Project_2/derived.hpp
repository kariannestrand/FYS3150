#ifndef DERIVED_CLASS_HPP
#define DERIVED_CLASS_HPP

#include "class.hpp" //Include the superclass header file here.

class DerivedClass : public MyClass { //Inherits all public and protected from SuperClass
private:
    /* data */
    double c0_, c1_, c2_;

public:
    DerivedClass (double c0, double c1, double c2);
    virtual ~DerivedClass (); //Destructor
    double compute_val(double x); //Computes the value y of the polynomial
};

#endif
