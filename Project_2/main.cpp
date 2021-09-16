// #include "class.hpp"
// #include "derived.hpp"

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace arma;
using namespace std;

mat anal(double N);
vec eigen_values(double N, double a, double d);

int main() {
    // Trying to do problem 3, but none of it is done with classes yet
    double n = 7;
    double N = (n-1);
    double h = 1/n;
    mat A = mat(N, N).fill(0.); // unfilled matrix
    double a = -1/(h*h);
    double d = 2/(h*h);
    vec lambda = vec(N);

    // filling main diagonal with 2/(h*h)
    for (int i = 0; i < N; i++){
            A(i, i) = d;
    }

    // filling sub- and superdiagonal with -1/(h*h)
    for (int i = 0; i < N-1; i++){
            A(i, i+1) = a;
            A(i+1, i) = a;
    }

    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, A);
    eigvec.print();
    cout << " " << endl;
    eigval.print();

    cout << " " << endl;
    //mat V = anal(N);
    mat V = normalise(anal(N), 2, 0);
    vec Lambda = eigen_values(N, a, d);
    V.print();
    cout << " " << endl;
    Lambda.print();




    /*
    // Line example
    double x = 2; //Point to evaluate the polynomials in
    double c0 = 1, c1 = 1; //Coefficients of the straight line.
    MyClass my_line = Line(c0, c1); //Create a Line object called my_line
    double y = my_line.compute_val(x); //Compute y = c0 + c1*x for a given x.
    */

    // quadratic example
    /*
    double x = 2; //Point to evaluate the polynomials in
    double c0 = 1, c1 = 1, c2 = 1;
    DerivedClass my_quad = DerivedClass(c0, c1, c2); //Call constructor and create Quadratic object my_quad.
    double y = my_quad.compute_val(x); //Compute y = c0 + c1*x + c2*x*x.
    */

    return 0;
}

mat anal(double N){
    mat V = mat(N, N).fill(0.);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            //V(j,i) = (i+1)*(j+1);
            V(j, i) = sin((i+1)*(j+1)*M_PI/(N+1));
        }
    }
    return V;
}

vec eigen_values(double N, double a, double d){
    vec Lambda = vec(N);
    for (int i = 0; i < N; i++){
        Lambda(i) = d + 2*a*cos((i+1)*M_PI/(N+1));
    }
    return Lambda;
}
