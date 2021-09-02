#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>

using namespace std; 

int main() {
    ofstream myfile; //trying to write the output to file
    myfile.open("problem2.txt");

   double x[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}; //array of xvalues
    for (int i=0; i < 11; i++) { //for loop to iterate over the function 
      myfile << "u(x)" << 1.0-(1.0-exp(-10.0))*x[i]-exp(-10.0*x[i]) << "x=" << x[i]; 
   }
   myfile.close();
}
