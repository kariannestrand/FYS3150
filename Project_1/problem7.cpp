#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;
     
// Declaring functions
void forward(vec, vec*, vec, vec*, int);
void backward(vec, vec, vec, vec*, int);
void elimination(vec, vec, vec, vec, vec*, int);
void writetofile(vec, vec, int);


int main(int argc, char* argv[]){
    int n = atoi(argv[1]);             // length of array
    double h = 1./n;                    // step size
    vec a = vec(n-1).fill(-1.);        // defining a-vector and filling with -1's
    vec b = vec(n).fill(2.);           // defining b-vector and filling with 2's
    vec c = vec(n-1).fill(-1.);        // defining c-vector and filling with -1's
    vec x = arma::linspace(0, 1, n);   // defining array x in [0, 1]
    vec f = 100*exp(-10*x);            // defining source term
    vec g = h*h*f;                     // defining solution vector g
    vec v = vec(n).fill(0.);           // creating empty solution vector v

    elimination(a, b, c, g, &v, n);    // calling elimination function
    writetofile(x, v, n);

    return 0;
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

// Writing output to file
void writetofile(vec x, vec v, int n){
  
    ofstream myfile;
    string filename = "problem7_" + to_string(n) + ".txt";
    myfile.open (filename);
    myfile << setw(15) << scientific << "x";
    myfile << setw(15) << scientific << "v" << endl;
    for (int i = 0; i < n; i++){
        myfile << setw(15) << scientific << x(i);
        myfile << setw(15) << scientific << v(i) << endl; // formatting txt-file
    }
    myfile.close();
    
}


