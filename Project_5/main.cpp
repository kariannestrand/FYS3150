#include "functions.hpp"

using namespace arma;
using namespace std;

int main(int argc, char const *argv[]){
    int M = 5;              // size of one side of outer matrix
    double h = 0.1;         // step size in x and y direction
    double dt = 0.1;        // step size for t
    int T = 100;            // total time

    cx_mat U_in = cx_mat(M-2, M-2, fill::randu);
    cx_vec u = state(U_in);


    cx_mat A = cx_mat((M-2)*(M-2), (M-2)*(M-2), fill::zeros);
    cx_vec a = cx_vec((M-2)*(M-2), fill::randu);
    double r = 2;


    A.submat(1, 1, M-2, M-2) = cx_mat(M-2, M-2, fill::ones);
    cout << A << endl;

    /*
    cout << a << endl;
    a.resize(M-2);
    cout << a << endl;
    */

    return 0;
}
