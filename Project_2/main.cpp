#include "class.hpp"

using namespace arma;
using namespace std;

//mat num(int N, double a, double d);
// mat eigen_vectors(int N);
<<<<<<< HEAD
vec eigen_values(int N, double a, double d);
double max_offdiag_symmetric(const mat &A, int &k, int &l);
void rotation(int N, mat &A, mat &R, int k, int l, double tol);
=======
// vec eigen_values(int N, double a, double d);
// double max_offdiag_symmetric(const mat &A, int &k, int &l);
// void rotation(int N, mat &A, mat &R, double k, double l, double tol);
>>>>>>> 92712d1d02e799a3297a79975999203836e47b4b

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


<<<<<<< HEAD
    //cout << " " << endl;
    //eigval.print();
    //eigvec.print()
    cout << " " << endl;
=======
    mat A = myclass.num();
>>>>>>> 92712d1d02e799a3297a79975999203836e47b4b

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
<<<<<<< HEAD
    //A.print();
    R.print();
    //cout << " " << endl;
    //normalise(R).print();
    cout << " " << endl;
    eigvec.print();
    //cout << " " << endl;
    //eigval.print();
    //cout << " " << endl;
    //A.print();
    //cout << count << endl;
    return 0;
}


double max_offdiag_symmetric(const mat &A, int &k, int &l){
    int n = A.n_rows;                       // declears length of matrix
    assert(A.is_square());                  // assert is an included function that makes sure the argument is true. We make sure the matrix is a nxn-matrix (square matrix)

    double maxval = 0.;                     // declears maxval
    for (int i = 0; i < n; i++){            // the for loop iterates in the lower triangle of the matrix where the row index always is always bigger than the column index
        for (int j = 0; j < n; j++){
            if (abs(A(i, j)) > maxval && i !=j){     // if the next element is bigger than the previous, save it as maxval
                maxval = abs(A(i, j));
                k = i;                      // indices of row for maxval
                l = j;                      // indices of column for maxval
            }
        }
    }
    return maxval;
}

/*
mat num(int N, double a, double d){
    mat A = mat(N, N).fill(0.);

    for (int i = 0; i < N-1; i++){
        // filling main diagonal with d:
        A(i, i) = d;
        // filling sub- and superdiagonal with a:
        A(i, i+1) = a;
        A(i+1, i) = a;
    }

    // filling the last element in the main diagonal with d:
    A(N-1, N-1) = d;

    return A;
}

mat eigen_vectors(int N){
    mat V = mat(N, N).fill(0.);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            //V(j,i) = (i+1)*(j+1);
            V(j, i) = sin((i+1)*(j+1)*M_PI/(N+1));
        }
    }
    return V;
}
*/

vec eigen_values(int N, double a, double d){
    vec Lambda = vec(N);
    for (int i = 0; i < N; i++){
        Lambda(i) = d + 2*a*cos((i+1)*M_PI/(N+1));
    }
    return Lambda;
}

void rotation(int N, mat &A, mat &R, int k, int l, double tol){
    double t,c,s,tau;

    if (A(k,l) != 0.0){
        tau = (A(l, l) - A(k, k))/(2*A(k, l));
        if (tau >= 0.0){
            t = 1/(tau + sqrt(1 + tau*tau));
        }else{
            t = -1/(-tau + sqrt(1 + tau*tau));
        }
        c = 1/sqrt(1 + t*t);
        s = c*t;
    }else{
        c = 1.0;
        s = 0.0;
    }

    double a_kk_m = A(k, k);
    A(k, k) = A(k, k)*c*c - 2*A(k, l)*c*s + A(l, l)*s*s;
    A(l, l) = A(l, l)*c*c + 2*A(k, l)*c*s + a_kk_m*s*s;
    A(k, l) = 0.0;
    A(l, k) = 0.0;


    // make sure to keep A_m(i, k) and A_m+1(i, k) separate! ?
    for (int i = 0; i < N; i++){
        if (i != l && i != k){
            double a_ik_m = A(i, k);
            A(i, k) = A(i, k)*c - A(i, l)*s;
            A(k, i) = A(i, k);
            A(i, l) = A(i, l)*c + a_ik_m*s;
            A(l, i) = A(i, l);
        }
    }

    for (int i = 0; i < N; i++){
        double r_ik_m = R(i, k);
        R(i, k) = R(i, k)*c - R(i, l)*s;
        R(i, l) = R(i, l)*c + r_ik_m*s;
    }
}
=======



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
>>>>>>> 92712d1d02e799a3297a79975999203836e47b4b
