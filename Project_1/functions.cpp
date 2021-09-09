#include "functions.hpp"


using namespace arma;
using namespace std;

// Calculating exact solution
vec exact(vec x){
  return 1 - (1 - exp(-10))*x - exp(-10*x);
}

// Forward substitution
void forward(vec a, vec *b, vec c, vec *g, int n){
    for (int i = 1; i < n-1; i++){
        (*b)(i) = (*b)(i)-a(i)/(*b)(i-1)*c(i-1);
        (*g)(i) = (*g)(i)-a(i)/(*b)(i-1)*(*g)(i-1);
    }
}

// Backwards substitution
void backward(vec b, vec c, vec g, vec *v, int n){
    (*v)(n-1) = g(n-1)/b(n-1);
    for (int i = n-2; i >= 0; i--){
      (*v)(i) = (g(i)-c(i)*(*v)(i+1))/b(i);
    }
}

// Gaussian elimination
void elimination(vec a, vec b, vec c, vec g, vec *v, int n){
    forward(a, &b, c, &g, n);
    backward(b, c, g, v, n);
}

// Writing to file exact.txt
void writetofile_exact(vec x, vec u){
    ofstream myfile;
    myfile.open ("exact.txt");
    myfile << setw(15) << scientific << "x";
    myfile << setw(15) << scientific << "u" << endl;
    for (int i = 0; i < 11; i++){
        myfile << setw(15) << scientific << x(i);
        myfile << setw(15) << scientific << u(i) << endl; // formatting txt-file
    }
    myfile.close();
}

// Writing to file approx_n.txt
void writetofile_approx(vec x, vec v, int n){
    ofstream myfile;
    string filename = "approx_" + to_string(n) + ".txt";
    myfile.open (filename);
    myfile << setw(15) << scientific << "x";
    myfile << setw(15) << scientific << "v" << endl;
    for (int i = 0; i < n; i++){
        myfile << setw(15) << scientific << x(i);
        myfile << setw(15) << scientific << v(i) << endl; // formatting txt-file
    }
    myfile.close();
    
}