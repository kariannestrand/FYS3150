#include "class.hpp"

using namespace arma;
using namespace std;

//This is how a constructor looks like in practice in the source file
MyClass::MyClass(int N, double a, double d){
  //definition of constructor
  //Assign member variables c0_ and c1_ to input variables c0 and c1, respectively.
  N_ = N;
  a_ = a;
  d_ = d;
}

//Definitions of other class methods come here.


mat MyClass::num(){
    mat A = mat(N_, N_).fill(0.);

    for (int i = 0; i < N_-1; i++){
        // filling main diagonal with d:
        A(i, i) = d_;
        // filling sub- and superdiagonal with a:
        A(i, i+1) = a_;
        A(i+1, i) = a_;
    }

    // filling the last element in the main diagonal with d:
    A(N_-1, N_-1) = d_;

    return A;
}

mat MyClass::eigen_vectors(){
    mat V = mat(N_, N_).fill(0.);
    for (int i = 0; i < N_; i++){
        for (int j = 0; j < N_; j++){
            //V(j,i) = (i+1)*(j+1);
            V(j, i) = sin((i+1)*(j+1)*M_PI/(N_+1));
        }
    }
    return V;
}
