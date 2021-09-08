#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;
double n = 10; // length of array
double h = 1/n; // step size

int main(){
  vec a = vec(n-1).fill(-1.); // defining a-vector
  vec b = vec(n).fill(2.); // defining b-vector 
  vec c = vec(n-1).fill(-1.); // defining c-vector
  vec x = arma::linspace(0, 1, n); // defining array x in [0, 1]
  vec f = 100*exp(-10*x); // defining source term
  vec g = h*h*f; // defining solution-vector g
  vec v = vec(n);
  for (int i = 1; i < n-1; i++){
      b(i) = b(i)-a(i)/b(i-1)*c(i-1);
      g(i) = g(i)-a(i)*b(i-1)*g(i-1);
  } 
  
 cout << b << endl;
 cout << g << endl;
 cout << v << endl;

  v(n) = g(n)/b(n);
  //cout << v(n);
  //for (int i = n-1; i >= 0; i--){
      //v(i) = (g(i)-c(i)*v(i+1))/b(i);
      //cout << v(i) << endl;
  //} 
  
  return 0;
}