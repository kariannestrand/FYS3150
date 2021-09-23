#include "class.hpp"

using namespace arma;
using namespace std;

int main(){
    int n = 7;
    int N = (n-1);
    double h = 1./n;
    double a = -1./(h*h);
    double d = 2./(h*h);
    double tol = 1e-12;
    int k;
    int l;

    MyClass myclass = MyClass(N, a, d);


    mat A = myclass.num();

    mat B = mat(4, 4).eye();
    B(3, 0) = 0.5;
    B(0, 3) = 0.5;
    B(1, 2) = -0.7;
    B(2, 1) = -0.7;

    mat R = mat(N, N).eye();


    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, A);

    mat V = normalise(myclass.eigen_vectors(), 2, 0);
    vec Lambda = myclass.eigen_values();


    double maxval_B = myclass.max_offdiag_symmetric(B, k, l);
    double maxval_A = myclass.max_offdiag_symmetric(A, k, l);


    int count = 0;
    while (maxval_A >= tol){
        count ++;
        myclass.rotation(A, R, k, l);
        maxval_A = myclass.max_offdiag_symmetric(A, k, l);
    }



    eigvec.print();
    cout << " " << endl;
    V.print();
    cout << " " << endl;
    R.print();
    cout << " " << endl;

    eigval.print();
    cout << " " << endl;
    Lambda.print();
    cout << " " << endl;
    A.print();


    //cout << count << endl;

    return 0;
}
