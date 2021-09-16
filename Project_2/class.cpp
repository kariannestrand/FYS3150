#include "class.hpp"

using namespace arma;
using namespace std;

//This is how a constructor looks like in practice in the source file
MyClass::MyClass(double c0, double c1){
  //definition of constructor
  //Assign member variables c0_ and c1_ to input variables c0 and c1, respectively.
  c0_ = c0;
  c1_ = c1;
}

//Definitions of other class methods come here.


double MyClass::compute_val(double x){
    //Returns the value of a straight line for an input x.
    return c0_ + c1_*x;
}
