// problem 2 of project 1
#include <iostream>
#include <armadillo>

double n = 10;

int main(){
  arma::vec x = arma::linspace(0, 1, n);
  arma::vec u = 1 - (1 - exp(-10))*x - exp(-10*x);
  std::cout << x << u << "\n";
  return 0;
}
