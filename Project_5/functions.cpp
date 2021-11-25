#include "functions.hpp"

using namespace arma;
using namespace std;


cx_vec state(cx_mat U_in){
    cx_vec u = U_in.as_col();
    return u;
}

cx_mat matrix(double r, cx_vec a, cx_mat A){
    /*
    cx_mat A = diagmat(a);
    cx_mat sub_mat_1 = A.submat(0, 0, M-2, M-2, fill::(ones));

    cx_mat sub_mat_1 = diagmat(1);
    for(int i = 0; i < (M-2); i++){
        for (int j = 0; j < (M-2); j++){
            if (i == j){
                mat1(i, j) = a;
            }
        }
    }
    */
    return A;
}
