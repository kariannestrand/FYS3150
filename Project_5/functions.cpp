#include "functions.hpp"

using namespace arma;
using namespace std;


cx_mat initial(double mean_x, double mean_y, double var_x, double var_y, int M){
    vec x = linspace(0, 1, M-2);  // x vector from 0 to 1 with M-3 steps
    vec y = linspace(0, 1, M-2);  // y vector from 0 to 1 with M-3 steps

    cx_vec distr_x = cx_vec(M-2);
    cx_vec distr_y = cx_vec(M-2);
    cx_mat U_in = cx_mat(M-2, M-2);

    cx_double i = cx_double(0.0, 1.0);
    double k_x = 1;
    double k_y = 1;
    for (int i = 0; i < M-2; i++){
        distr_x(i) = 1/(var_x*sqrt(2*M_PI))*exp(-(x(i) - mean_x)/(var_x)*(x(i) - mean_x)/(var_x))*exp(i*k_x*x(i));
        distr_y(i) = 1/(var_y*sqrt(2*M_PI))*exp(-(y(i) - mean_y)/(var_y)*(y(i) - mean_y)/(var_y))*exp(i*k_y*y(i));
    }

    for (int i = 0; i < M-2; i++){
        for (int j = 0; j < M-2; j++){
            U_in(i, j) = distr_x(i)*distr_y(j);
        }
    }

    U_in = U_in/norm(U_in);  // making sure probability distribution starts out normalized to 1

    return U_in;
}


cx_mat potential(double v0, int M, int size_slit, int size_between_slit){
    cx_mat V = cx_mat(M-2, M-2, fill::zeros);
    for (int i = 0; i < M-2; i++){
        for (int j = 0; j < M-2; j++){
            if (j == (M-2)/2){
                V(i, j) = v0;

                for (int k = size_between_slit; k < size_between_slit + size_slit; k++){
                    if (i == (M-2)/2 + k){
                        V(i, j) = 0.0;
                    }
                    else if (i == (M-2)/2 - k){
                        V(i, j) = 0.0;
                    }
                }
            }
        }
    }
    return V;
}


void vector_ab(cx_double r, double dt, int M, cx_vec &a, cx_vec &b, cx_mat V){
    cx_vec v = V.as_col();
    complex<double> i(0.0, 1.0);

    for (int k = 0; k < pow(M-2, 2); k++){
        a(k) = 1. + 4.*r + i*dt/2.0 * v(k);
        b(k) = 1. - 4.*r - i*dt/2.0 * v(k);
    
    }

}


void matrix(cx_double r, cx_vec a, cx_vec b, sp_cx_mat &A, cx_mat &B, int M){
    for (int i = 0; i < (M-2)*(M-2); i++){
        A(i, i) = a(i);
        B(i, i) = b(i);
    }
    for (int i = 0; i < (M-2)*(M-2)-(M-2); i++){
        A(i, i + (M-2)) = -r;
        A(i + (M-2), i) = -r;
        B(i, i + (M-2)) = r;
        B(i + (M-2), i) = r;
    }
    for (int i = 0; i < (M-2)*(M-2)-1; i++){
        A(i, i + 1) = -r;
        A(i + 1, i) = -r;
        B(i, i + 1) = r;
        B(i + 1, i) = r;
        for (int j = 1; j < (M-2)*(M-2)-1; j++){
            if (i == j*(M-2) - 1){
                A(i, i + 1) = 0;
                A(i + 1, i) = 0;
                B(i, i + 1) = 0;
                B(i + 1, i) = 0;
            }
        }
    }
}


cx_vec CrankNicolson(cx_mat U_in, cx_mat B, sp_cx_mat A, int N_t){
    cx_vec u = U_in.as_col();
    cx_vec b = cx_vec(u.size());
    complex<double> p;

    
    for (int i = 0; i < N_t; i++){
        b = B*u;
        u = spsolve(A, b);
        p = cdot(u, u);
        cout << p << endl;
    }
    

    return u;
}
