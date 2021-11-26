#include "functions.hpp"

using namespace arma;
using namespace std;

int main(int argc, char const *argv[]){
    int M = 5;               // size of one side of outer matrix
    double h = 1./(M-1);           // step size in x and y direction
    double dt = 0.01;          // step size for t
    int N_t = 100;
    cx_double r = cx_double(0.0, dt/(2*h*h));
    double v0 = 10000;
    int size_slit = 1;
    int size_between_slit = 1;

    cx_vec a = cx_vec((M-2)*(M-2));
    cx_vec b = cx_vec((M-2)*(M-2));
    sp_cx_mat A = sp_cx_mat((M-2)*(M-2), (M-2)*(M-2));
    cx_mat B = cx_mat((M-2)*(M-2), (M-2)*(M-2), fill::zeros);

    double mean_x = 0.5;
    double mean_y = 0.5;
    double var_x = 0.25;
    double var_y = 0.25;

    cx_mat U_in = initial(mean_x, mean_y, var_x, var_y, M);      // creates u(0, 0, 0) matrix
    cx_mat V = potential(v0, M, size_slit, size_between_slit);
    vector_ab(r, dt, M, a, b, V);                                // creating a and b vector
    matrix(r, a, b, A, B, M);                                    // creating A and B matrix


    cx_vec U_in_vec = CrankNicolson(U_in, B, A, N_t);              // solving u(0, 0, 1)

    return 0;
}