#include "functions.hpp"
#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;


int main(int argc, char* argv[]){
    int n = atoi(argv[1]);                    // length of array
    double h = 1./n;                          // step size
    vec a = vec(n-1).fill(-1.);               // defining a-vector filled with -1's
    vec b = vec(n).fill(2.);                  // defining b-vector filled with 2's
    vec c = vec(n-1).fill(-1.);               // defining c-vector filled with -1's
    vec x = arma::linspace(0+h, 1-h, n);      // defining array x in [0, 1]
    vec f = 100*exp(-10*x);                   // defining source term
    vec g = h*h*f;                            // defining solution vector g
    vec v = vec(n).fill(0.);                  // creating empty solution vector v


    vec u = exact(x);                         // calling exact solution function
    elimination(a, b, c, g, &v, n);           // calling gaussian elimination function
    vec Delta = abs_err(u, v);                // calling absolute error function
    vec epsilon = rel_err(u, v);              // calling relative error function
    writetofile(x, u, v, Delta, epsilon, n);  // calling writing to file function

    return 0;
}
