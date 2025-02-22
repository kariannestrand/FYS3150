#include "jacobi_rotation.hpp"

using namespace arma;
using namespace std;

JacobiRotation::JacobiRotation(int N, double a, double d){
  N_ = N;
  a_ = a;
  d_ = d;
}

mat JacobiRotation::tridiag_matrix(){
    mat A = mat(N_, N_).fill(0.);      // declearing NxN-matrix filled with zeros

    for (int i = 0; i < N_-1; i++){
        // filling main diagonal with d:
        A(i, i) = d_;
        // filling sub- and superdiagonal with a:
        A(i, i+1) = a_;
        A(i+1, i) = a_;
    }

    // filling the last element in the main diagonal with d:
    A(N_-1, N_-1) = d_;

    return A;
}

mat JacobiRotation::eigen_vectors(){
    mat V = mat(N_, N_).fill(0.);      // declearing NxN-matrix filled with zeros
    for (int i = 0; i < N_; i++){
        for (int j = 0; j < N_; j++){
            V(j, i) = sin((i+1)*(j+1)*M_PI/(N_+1));    // filling matrix with analytical expression for eigenvectors
        }
    }
    return V;
}

vec JacobiRotation::eigen_values(){
    vec Lambda = vec(N_);             // declearing NxN-matrix filled with zeros
    for (int i = 0; i < N_; i++){
        Lambda(i) = d_ + 2*a_*cos((i+1)*M_PI/(N_+1));  // filling matrix with analytical expression for eigenvalues
    }
    return Lambda;
}

double JacobiRotation::max_offdiag_symmetric(mat &A, int &k, int &l){
    int n = A.n_rows;                             // declears length of matrix
    assert(A.is_square());                        // assert is an included function that makes sure the argument is true. We make sure the matrix is a nxn-matrix (square matrix)

    double maxval = 0.;                            // declears maxval
    for (int j = 0; j < n; j++){                   // the for loop iterates in the lower triangle of the matrix where the row index always is bigger than the column index
        for (int i = j + 1; i < n; i++){
            if (abs(A(i, j)) > maxval && i !=j){  // if the next element is bigger than the previous, save it as maxval
                maxval = abs(A(i, j));
                k = i;                            // indices of row for maxval
                l = j;                            // indices of column for maxval
            }
        }
    }
    return maxval;
}

void JacobiRotation::rotation(mat &A, mat &R, int k, int l){
    // Jacobi´s rotation method:
    double t, c, s, tau;

    if (A(k,l) != 0.0){
        tau = (A(l, l) - A(k, k))/(2*A(k, l));
        if (tau >= 0.0){
            t = 1/(tau + sqrt(1 + tau*tau));
        }else{
            t = -1/(-tau + sqrt(1 + tau*tau));
        }
        c = 1/sqrt(1 + t*t);
        s = c*t;
    }
    else{
        c = 1.0;
        s = 0.0;
    }

    double a_kk_m = A(k, k);         // makes sure to keep A^m_(i, k) and A^(m+1)_(i, k) separate
    A(k, k) = A(k, k)*c*c - 2*A(k, l)*c*s + A(l, l)*s*s;
    A(l, l) = A(l, l)*c*c + 2*A(k, l)*c*s + a_kk_m*s*s;
    A(k, l) = 0.0;
    A(l, k) = 0.0;


    for (int i = 0; i < N_; i++){
        if (i != l && i != k){
            double a_ik_m = A(i, k); // makes sure to keep A^m_(i, k) and A^(m+1)_(i, k) separate
            A(i, k) = A(i, k)*c - A(i, l)*s;
            A(k, i) = A(i, k);
            A(i, l) = A(i, l)*c + a_ik_m*s;
            A(l, i) = A(i, l);
        }
    }

    for (int i = 0; i < N_; i++){
        double r_ik_m = R(i, k);    // makes sure to keep R^m_(i, k) and R^(m+1)_(i, k) separate
        R(i, k) = R(i, k)*c - R(i, l)*s;
        R(i, l) = R(i, l)*c + r_ik_m*s;
    }
}

void JacobiRotation::write(mat R, int number){
    // function that writes eigenvectors found with Jacobi´s rotation method to file
    ofstream file;
    file.open("eigvec_" + to_string(number) + ".txt");
    for(int i = 0; i < R.size(); ++i){
        file << R[i] << endl;
    }
    file.close();
}
