#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std; 

int main() {
    //array of xvalues
   double ar[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    for (int i=0; i < 11; i++) {
      cout << "u=" << 1.0-(1.0-exp(-10.0))*ar[i]-exp(-10.0*ar[i]);
   }
}
