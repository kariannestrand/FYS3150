#include "functions.hpp"

using namespace arma;
using namespace std;

int main(int argc, char const *argv[]){
    bool write = true;                             // writes internal state cube U(x, y, t) to bin-file in true
    string name = "T_0002.bin";                    // name of bin-file

    double h = 0.005;                              // step size in x and y direction
    double M = 1.0/h + 1.0;                        // size of one side of outer matrix
    double dt = 2.5*pow(10, -5);                   // step size for t
    double T = 0.002;                              // total time
    double N = T/dt;                               // number of time steps
    double v0 = pow(10, 10);                       // constant potential value inside barriers
    cx_double r = cx_double(0.0, dt/(2*h*h));
    int size_slit = 1;
    int size_between_slit = 1;


    cx_vec a = cx_vec((M-2)*(M-2));                                      // initializing a vector
    cx_vec b = cx_vec((M-2)*(M-2));                                      // initializing b vector
    sp_cx_mat A = sp_cx_mat((M-2)*(M-2), (M-2)*(M-2));                   // initializing A matrix
    cx_mat B = cx_mat((M-2)*(M-2), (M-2)*(M-2), fill::zeros);            // initializing B matrix

    double mean_x = 0.25;                                                // x_c
    double mean_y = 0.5;                                                 // y_c
    double var_x = 0.05;                                                 // sigma_x
    double var_y = 0.2;                                                  // sigma_y
    double p_x = 200.0;                                                  // p_x
    double p_y = 0.0;                                                    // p_y

    cx_mat U = initial(mean_x, mean_y, var_x, var_y, p_x, p_y, M, h);    // initializing internal state matrix U(x, y, 0)
    cx_mat V = potential(v0, M, size_slit, size_between_slit);           // generating potential matrix V(x, y)
    vector_ab(r, dt, M, a, b, V);                                        // generating a and b vector
    matrix(r, a, b, A, B, M);                                            // generating A and B matrix
    cx_cube U_cube = CrankNicolson(U, B, A, M, N);                       // generating internal state cube U(x, y, t) where each slize is a time step

    if (write){
        write_to_file(U_cube, name);                                     // writes U_cube to bin-file
    }

    return 0;
}
