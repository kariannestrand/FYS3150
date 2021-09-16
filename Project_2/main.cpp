#include "class.hpp"
#include "derived.hpp"

using namespace arma;
using namespace std;

int main() {
    // Line example
    double x = 2; //Point to evaluate the polynomials in
    double c0 = 1, c1 = 1; //Coefficients of the straight line.
    MyClass my_line = Line(c0, c1); //Create a Line object called my_line
    double y = my_line.compute_val(x); //Compute y = c0 + c1*x for a given x.

    // quadratic example
    /*
    double x = 2; //Point to evaluate the polynomials in
    double c0 = 1, c1 = 1, c2 = 1;
    DerivedClass my_quad = DerivedClass(c0, c1, c2); //Call constructor and create Quadratic object my_quad.
    double y = my_quad.compute_val(x); //Compute y = c0 + c1*x + c2*x*x.
    */

    return 0;
}
