#include "functions.hpp"

using namespace arma;
using namespace std;

int main(int argc, char const *argv[]){
    int M = 6;              // size of one side of outer matrix
    double h = 0.1;         // step size in x and y direction
    double dt = 0.1;        // step size for t
    int T = 100;            // total time
    double r = 2;

    cx_mat U_in = cx_mat(M-2, M-2, fill::randu);

    cx_vec a = cx_vec((M-2)*(M-2));
    cx_vec b = cx_vec((M-2)*(M-2));

    cx_mat A = cx_mat((M-2)*(M-2), (M-2)*(M-2), fill::zeros);
    cx_mat B = cx_mat((M-2)*(M-2), (M-2)*(M-2), fill::zeros);

    vector_ab(r, dt, M, a, b);
    matrix(r, a, b, A, B, M);

    solver(U_in, B, A);

    return 0;
}