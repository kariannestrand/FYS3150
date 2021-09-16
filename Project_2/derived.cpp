#include "derived.hpp"

DerivedClass::DerivedClass(double c0, double c1, double c2) : Line(c0, c1){ //Here we reuse the constructor of the Line class. Note the syntax!
    c2_ = c2; //Assign the higher order coefficient. The others are assigned in Line(c0,c1).
}

double DerivedClass::compute_val(double x){
    //Computes the quadratic polynomial value y at x.
    //Reuses Line::compute_val(x) to compute the contribution from the straight line and tacks on the quadratic explicitly.
    return MyClass::compute_val(x) + c2_*x*x;
}
