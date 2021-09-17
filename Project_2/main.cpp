// #include "class.hpp"
// #include "derived.hpp"

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace arma;
using namespace std;

mat num(double N, double a, double d);
mat anal(double N);
vec eigen_values(double N, double a, double d);
double max_offdiag_symmetric(const arma::mat& B, int& k, int& l);

int main() {
    double n = 7;
    double N = (n-1);
    double h = 1/n;
    double a = -1/(h*h);
    double d = 2/(h*h);

    mat A = num(N, a, d);
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, A);

    mat V = normalise(anal(N), 2, 0);
    vec Lambda = eigen_values(N, a, d);

    eigvec.print();
    cout << " " << endl;
    eigval.print();
    cout << " " << endl;
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


// A function that finds the max off-diag element of a symmetric matrix A.
// - The matrix indices of the max element are returned by writing to the
//   int references k and l (row and column, respectively)
// - The value of the max element A(k,l) is returned as the function
//   return value
double max_offdiag_symmetric(const arma::mat& B, int& k, int& l){
    mat B = mat(N, N).fill(0.);
    for i in range
  // Get size of the matrix A. Use e.g. A.n_rows, see the Armadillo documentation

  // Possible consistency checks:
  // Check that A is square and larger than 1x1. Here you can for instance use A.is_square(),
  // see the Armadillo documentation.
    cout << B.is_square() << endl;
  //
  // The standard function 'assert' from <assert.h> can be useful for quick checks like this
  // during the code development phase. Use it like this: assert(some condition),
  // e.g assert(a==b). If the condition evaluates to false, the program is killed with
  // an assertion error. More info: https://www.cplusplus.com/reference/cassert/assert/

  // Initialize references k and l to the first off-diagonal element of A

  // Initialize a double variable 'maxval' to A(k,l). We'll use this variable
  // to keep track of the largest off-diag element.

  // Loop through all elements in the upper triangle of A (not including the diagonal)
  // When encountering a matrix element with larger absolute value than the current value of maxval,
  // update k, l and max accordingly.

  return maxval; // Return maxval 
}

mat num(double N, double a, double d){
    mat A = mat(N, N).fill(0.);
    for (int i = 0; i < N; i++){
        // filling main diagonal with d:
        A(i, i) = d;
    }

    for (int i = 0; i < N-1; i++){
        // filling sub- and superdiagonal with a:
        A(i, i+1) = a;
        A(i+1, i) = a;
    }
    return A;
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
