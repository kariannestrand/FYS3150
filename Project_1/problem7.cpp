#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;


vec forward(vec, vec, vec, vec, int);
vec backward(vec, vec, vec, vec, int);
vec elimination(vec, vec, vec, vec, vec, int);
//void writetofile(vec, vec, int);


int main(int argc, char* argv[]){
    int n = atoi(argv[1]);             // length of array
    double h = 1/n;                    // step size
    vec a = vec(n-1).fill(-1.);        // defining a-vector and filling with -1's
    vec b = vec(n).fill(2.);           // defining b-vector and filling with 2's
    vec c = vec(n-1).fill(-1.);        // defining c-vector and filling with -1's
    vec x = arma::linspace(0, 1, n);   // defining array x in [0, 1]
    vec f = 100*exp(-10*x);            // defining source term
    vec g = h*h*f;                     // defining solution vector g
    vec v = vec(n);                    // creating empty solution vector v

    v = elimination(a, b, c, g, v, n);

/*
    ofstream myfile;
    myfile.open ("problem7_10000.txt");
    myfile << setw(15) << setprecision(8) << "x";
    myfile << setw(15) << setprecision(8) << "v" << endl;
    for (int i = 0; i < n; i++){
        myfile << setw(15) << setprecision(8) << x(i);
        myfile << setw(15) << setprecision(8) << v(i) << endl; // formatting txt-file
    }
    myfile.close();
*/
    return 0;
}


vec forward(vec a, vec b, vec c, vec g, int n){
    //Forward substitution
    for (int i = 1; i < n-1; i++){
        b(i) = b(i)-a(i)/b(i-1)*c(i-1);
        g(i) = g(i)-a(i)/b(i-1)*g(i-1);
    }
    return b, g;
}
vec backward(vec b, vec c, vec g, vec v, int n){
    //Backwards substitution
    v(n-1) = g(n-1)/b(n-1);
    for (int i = n-2; i >= 0; i--){
      v(i) = (g(i)-c(i)*v(i+1))/b(i);
    }
    return v;
}
vec elimination(vec a, vec b, vec c, vec g, vec v, int n){
    b, g = forward(a, b, c, g, n);
    v = backward(b, c, g, v, n);
    return v;
}
/*
void writetofile(vec x, vec v, int n){

}
*/
