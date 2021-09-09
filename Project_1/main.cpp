#include "functions.hpp"
#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;
     

int main(int argc, char* argv[]){
    int n = atoi(argv[1]);             // length of array
    double h = 1./n;  
    vec x_exact = arma::linspace(0, 1, 11);                  // step size
    vec a = vec(n-1).fill(-1.);        // defining a-vector and filling with -1's
    vec b = vec(n).fill(2.);           // defining b-vector and filling with 2's
    vec c = vec(n-1).fill(-1.);        // defining c-vector and filling with -1's
    vec x_approx = arma::linspace(0+h, 1-h, n);   // defining array x in [0, 1]
    vec f = 100*exp(-10*x_approx);            // defining source term
    vec g = h*h*f;                     // defining solution vector g
    vec v = vec(n).fill(0.);           // creating empty solution vector v


    vec u = exact(x_exact);
    elimination(a, b, c, g, &v, n);    // calling elimination function
    writetofile_exact(x_exact, u);
    writetofile_approx(x_approx, v, n);

    return 0;
}

