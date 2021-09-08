#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;
double n = 10;       // length of array
double h = 1/n;      // step size

int main(){
<<<<<<< HEAD
  vec a = vec(n-1).fill(-1.);          // defining a-vector
  vec b = vec(n).fill(2.);             // defining b-vector
  vec c = vec(n-1).fill(-1.);          // defining c-vector
  vec x = arma::linspace(0, 1, n);     // defining array x in [0, 1]
  vec f = 100*exp(-10*x);              // defining source term
  vec g = h*h*f;                       // defining solution-vector g
  vec v = vec(n);
  for (int i = 1; i < n-1; i++){
      b(i) = b(i)-a(i)/b(i-1)*c(i-1);
      g(i) = g(i)-a(i)*b(i-1)*g(i-1);
  }

  v(n-1) = g(n-1)/b(n-1);
  cout << v << endl;
  //for (int i = n-2; i >= 0; i--){
  //    v(i) = (g(i)-c(i)*v(i+1))/b(i);
  //    cout << v(i) << endl;
  //}
=======
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

>>>>>>> dd70c7191804af4aa788c8b3bcdee1e8c17bb653

    //cout << g << endl;
    //cout << b << endl;
    v(n-1) = g(n-1)/b(n-1);
    //cout << v << endl;
    for (int i = n-2; i >= 0; i--){
      v(i) = (g(i)-c(i)*v(i+1))/b(i);
    }
    cout << v << endl;
    //cout << b << endl;
    //cout << g << endl;
    return 0;
}
