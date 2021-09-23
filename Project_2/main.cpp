#include "class.hpp"
#include <vector>

using namespace arma;
using namespace std;

int main(int argc, char* argv[]){
    int n = atoi(argv[1]);              // declearing number of points
    int N = (n-1);                      // declearing number of steps
    double h = 1./n;                    // declearing step size
    double a = -1./(h*h);               // declearing sub- and superdiagonal values
    double d = 2./(h*h);                // declearing main diagonal values
    double tol = 1e-12;                 // declearing tolerance
    int k;                              // declearing maximum value row-index
    int l;                              // declearing maximum value column-index

    MyClass myclass = MyClass(N, a, d); // calling MyClass

    mat A = myclass.tridiag_matrix();   // creating matrix A with a and d on tridiagonal

    mat B = mat(4, 4).eye();            // creating 4x4-matrix B
    B(3, 0) = 0.5;
    B(0, 3) = 0.5;
    B(1, 2) = -0.7;
    B(2, 1) = -0.7;

    mat R = mat(N, N).eye();


    // finding eigenvectors and eigenvalues of A with numerical method
    vec eigval;                         // declearing eigenvalues
    mat eigvec;                         // declearing eigenvectors
    eig_sym(eigval, eigvec, A);

    // finding eigenvectors and eigenvalues of A with analytical method
    mat V = normalise(myclass.eigen_vectors(), 2, 0);     // eigenvector
    vec Lambda = myclass.eigen_values();                  // eigenvalues


    double maxval_B = myclass.max_offdiag_symmetric(B, k, l);
    double maxval_A = myclass.max_offdiag_symmetric(A, k, l);


    int count = 0;
    while (maxval_A >= tol){
        count ++;
        myclass.rotation(A, R, k, l);
        maxval_A = myclass.max_offdiag_symmetric(A, k, l);
    }


    // prints eigenvalues and eigenvectors of A for different methods if eig = true
    bool print_eig = false;
    if (print_eig){
        cout << " " << endl;
        cout << "Numerical eigenvectors of A:" << endl;
        eigvec.print();
        cout << endl;
        cout << "Analytical eigenvectors of A:" << endl;
        V.print();
        cout << " " << endl;
        cout << "Eigenvectors of A with Jacobi´s rotation method:" << endl;
        R.print();
        cout << " " << endl;

        cout << "Numerical eigenvalues of A:" << endl;
        eigval.print();
        cout << " " << endl;
        cout << "Analytical eigenvalues of A:" << endl;
        Lambda.print();
        cout << " " << endl;
        cout << "Eigenvalues of A with Jacobi´s rotation method:" << endl;
        A.print();
    }


    // prints number of required transformations if print_count = true
    bool print_count = false;
    if (print_count){
        cout << count << endl;
    }

    // makes txt-files with R-columns if make_txt_files = true
    bool make_txt_files = true;
    if (make_txt_files){
        vec A_eigvals = vec(N).fill(0.);
        for (int i = 0; i < N; i++){
            A_eigvals(i) = A(i, i);
        }
        uvec eigval_indices = sort_index(A_eigvals);

        myclass.write(R.col(eigval_indices(0)), 1);
        myclass.write(R.col(eigval_indices(1)), 2);
        myclass.write(R.col(eigval_indices(2)), 3);
    }
    return 0;
}
