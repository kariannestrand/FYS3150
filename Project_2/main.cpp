// #include "class.hpp"
// #include "derived.hpp"

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <assert.h>

using namespace arma;
using namespace std;

mat num(double N, double a, double d);
mat anal(double N);
vec eigen_values(double N, double a, double d);
double max_offdiag_symmetric(const mat& B, int& k, int& l);

int main() {
    int n = 7;
    int N = (n-1);
    double h = 1./n;
    double a = -1./(h*h);
    double d = 2./(h*h);

    mat A = num(N, a, d);
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, A);

    mat V = normalise(anal(N), 2, 0);
    vec Lambda = eigen_values(N, a, d);

    /*
    eigvec.print();
    cout << " " << endl;
    eigval.print();
    cout << " " << endl;
    V.print();
    cout << " " << endl;
    Lambda.print();
    */

    int s = 4;
    mat B = mat(s, s).fill(0.);
    B.eye();
    B(3, 0) = 0.5;
    B(0, 3) = 0.5;
    B(1, 2) = -0.7;
    B(2, 1) = -0.7;

    int k;
    int l;

    double maxval = max_offdiag_symmetric(B, k, l);
    cout << k <<  " " << l << endl;
    cout << maxval << endl;




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


double max_offdiag_symmetric(const mat &B, int &k, int &l){
    int n = B.n_rows;                       // declears length of matrix
    assert(B.is_square());                  // assert is an included function that makes sure the argument is true. We make sure the matrix is a nxn-matrix (square matrix)

    double maxval = 0.;                     // declears maxval
    for (int j = 0; j < n; j++){            // the for loop iterates in the lower triangle of the matrix where the row index always is always bigger than the column index
        for (int i = j+1; i < n; i++){
            if (abs(B(i, j)) > maxval){     // if the next element is bigger than the previous, save it as maxval
                maxval = abs(B(i, j));
                k = i;                      // indices of row for maxval
                l = j;                      // indices of column for maxval
            }
        }
    }
    return maxval;
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
