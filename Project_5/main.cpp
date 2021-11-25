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

    /*
    cx_mat A = cx_mat((M-2)*(M-2), (M-2)*(M-2), fill::zeros);
    cx_vec a = cx_vec((M-2)*(M-2), fill::randu);
    */
    double r = 2;


<<<<<<< HEAD
    A.submat(0, 0, M-3, M-3) = cx_mat(M-2, M-2, fill::ones);
=======
    mat A = mat((M-2)*(M-2), (M-2)*(M-2), fill::zeros);

    for (int i = 0; i < (M-3)*(M-3); i++){
        // filling main diagonal with a-vector:
        //A(i, i) = a(i);
        A(i,i) = 1.0;
        // filling sub- and superdiagonal with a:
        A(i, i + (M-2)) = -r;
        A(i + (M-2), i) = -r;
    }

    //A = matdiag(a);


    //submat_1 = A.submat(0, 0, M-3, M-3) = cx_mat(M-2, M-2, fill::eye);
>>>>>>> 41ead1596e88414b32c870240e3daf599bdcd6f0
    cout << A << endl;

    /*
    cout << a << endl;
    a.resize(M-2);
    cout << a << endl;
    */

    return 0;
}
