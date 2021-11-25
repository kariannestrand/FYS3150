#include "functions.hpp"

using namespace arma;
using namespace std;


cx_vec state(cx_mat U_in){
    cx_vec u = U_in.as_col();
    return u;
}

void matrix(double r, cx_vec a, cx_vec b, cx_mat &A, cx_mat &B, int M){
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
