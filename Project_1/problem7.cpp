#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;

double n = 10000; // length of array
double h = 1/n; // step size

//defining forward function
//defining backwards function
//defining function for writing to file

int main(){
    vec a = vec(n-1).fill(-1.); // defining a-vector and filling with -1's
    vec b = vec(n).fill(2.); // defining b-vector and filling with 2's
    vec c = vec(n-1).fill(-1.); // defining c-vector and filling with -1's
    vec x = arma::linspace(0, 1, n); // defining array x in [0, 1]
    vec f = 100*exp(-10*x); // defining source term
    vec g = h*h*f; // defining solution vector g
    vec v = vec(n); // creating empty solution vector v
    //Forward substitution
    for (int i = 1; i < n-1; i++){
        b(i) = b(i)-a(i)/b(i-1)*c(i-1);
        g(i) = g(i)-a(i)/b(i-1)*g(i-1);
    }

    v(n-1) = g(n-1)/b(n-1);
    //Backwards substitution
    for (int i = n-2; i >= 0; i--){
      v(i) = (g(i)-c(i)*v(i+1))/b(i);
      //cout << v(i) << endl;
    }

    ofstream myfile;
    myfile.open ("problem7_10000.txt");
    myfile << setw(15) << setprecision(8) << "x";
    myfile << setw(15) << setprecision(8) << "v" << endl;
    for (int i = 0; i < n; i++){
        myfile << setw(15) << setprecision(8) << x(i);
        myfile << setw(15) << setprecision(8) << v(i) << endl; // formatting txt-file
    }
    myfile.close();
    return 0;
}
