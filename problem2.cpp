// problem 2 of project 1
#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;
double n = 11; // length of array

int main(){
  vec x = arma::linspace(0, 1, n); // defining array x in [0, 1]
  vec u = 1 - (1 - exp(-10))*x - exp(-10*x); // calculating solution

  // writing solution to file called problem2.txt
  ofstream myfile;
  myfile.open ("problem2.txt");
  myfile << setw(15) << setprecision(8) << "x";
  myfile << setw(15) << setprecision(8) << "u(x)" << endl;
  for (int i = 0; i < n; i++){
    myfile << setw(15) << setprecision(8) << x(i);
    myfile << setw(15) << setprecision(8) << u(i) << endl; // formatting txt-file
  }
  myfile.close();

  return 0;
}
