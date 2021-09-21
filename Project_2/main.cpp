// #include "class.hpp"
// #include "derived.hpp"

#include <armadillo>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <assert.h>

using namespace arma;
using namespace std;

mat num(int N, double a, double d);
mat eigen_vectors(int N);
vec eigen_values(int N, double a, double d);
double max_offdiag_symmetric(const mat &A, int &k, int &l);
void rotation(int N, mat &A, mat &R, double k, double l, double tol);

int main(){
    int n = 7;
    int N = (n-1);
    double h = 1./n;
    double a = -1./(h*h);
    double d = 2./(h*h);
    double tol = 1e-8;

    mat A = num(N, a, d);
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, A);

    // mat V = normalise(eigen_vectors(N), 2, 0);
    // vec Lambda = eigen_values(N, a, d);

    /*
    eigvec.print();
    cout << " " << endl;
    eigval.print();
    cout << " " << endl;
    V.print();
    cout << " " << endl;
    Lambda.print();
    */

    mat B = mat(4, 4).eye();
    B(3, 0) = 0.5;
    B(0, 3) = 0.5;
    B(1, 2) = -0.7;
    B(2, 1) = -0.7;

    int k;
    int l;

    double maxval_B = max_offdiag_symmetric(B, k, l);


    double maxval_A = max_offdiag_symmetric(A, k, l);
    cout << k <<  " " << l << endl;
    cout << maxval_A << endl;

    mat R = mat(N, N).eye();

    rotation(N, A, R, k, l, tol);
    A.print();
    R.print();

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


double max_offdiag_symmetric(const mat &A, int &k, int &l){
    int n = A.n_rows;                       // declears length of matrix
    assert(A.is_square());                  // assert is an included function that makes sure the argument is true. We make sure the matrix is a nxn-matrix (square matrix)

    double maxval = 0.;                     // declears maxval
    for (int j = 0; j < n; j++){            // the for loop iterates in the lower triangle of the matrix where the row index always is always bigger than the column index
        for (int i = j+1; i < n; i++){
            if (abs(A(i, j)) > maxval){     // if the next element is bigger than the previous, save it as maxval
                maxval = abs(A(i, j));
                k = i;                      // indices of row for maxval
                l = j;                      // indices of column for maxval
            }
        }
    }
    return maxval;
}

mat num(int N, double a, double d){
    mat A = mat(N, N).fill(0.);

    for (int i = 0; i < N-1; i++){
        // filling main diagonal with d:
        A(i, i) = d;
        // filling sub- and superdiagonal with a:
        A(i, i+1) = a;
        A(i+1, i) = a;
    }

    // filling the last element in the main diagonal with d:
    A(N-1, N-1) = d;

    return A;
}

mat eigen_vectors(int N){
    mat V = mat(N, N).fill(0.);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            //V(j,i) = (i+1)*(j+1);
            V(j, i) = sin((i+1)*(j+1)*M_PI/(N+1));
        }
    }
    return V;
}

vec eigen_values(int N, double a, double d){
    vec Lambda = vec(N);
    for (int i = 0; i < N; i++){
        Lambda(i) = d + 2*a*cos((i+1)*M_PI/(N+1));
    }
    return Lambda;
}

void rotation(int N, mat &A, mat &R, double k, double l, double tol){
    while (A(k, l) > tol){

        double tau = (A(l, l) - A(k, k))/(2*A(k, l));
        double t;
        double c = 1/sqrt(1 + tau*tau);
        double s = c*t;

        if (tau > 0){
            t = 1/(tau + sqrt(1 + tau*tau));
        }
        if (tau < 0){
            t = -1/(-tau + sqrt(1 + tau*tau));
        }
        
        /*
        if (A(k, l) == 0){
            c = 1;
            s = 0;
            t = 0;
        }
        */

        A(k, k) = A(k, k)*c*c - 2*A(k, l)*c*s + A(l, l)*s*s;
        A(l, l) = A(l, l)*c*c + 2*A(k, l)*c*s + A(k, k)*s*s;
        A(l, k) = 0;
        A(l, k) = 0;


        // make sure to keep A_m(i, k) and A_m+1(i, k) separate! ?
        for (int i = 0; i != l && i != k && i < N; i++){
            double a_ik_m = A(i, k);
            A(i, k) = A(i, l)*c - A(i, l)*s;
            A(k, i) = A(i, k);
            A(i, l) = A(i, l)*c + a_ik_m*s;
            A(l, i) = A(i, l);
        }

        for (int i = 0; i < N; i++){
            R(i, k) = R(i, k)*c - R(i, l)*s;
            R(i, l) = R(i, l)*c + R(i, k)*s;
        }
    }
}
