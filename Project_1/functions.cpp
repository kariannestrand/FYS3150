#include "functions.hpp"

using namespace arma;
using namespace std;



// Calculating exact solution
vec exact(vec x){
  return 1 - (1 - exp(-10))*x - exp(-10*x);
}


// Forward substitution, general algorithm
void forward_gen(vec a, vec *b, vec c, vec *g, int n){
    for (int i = 1; i < n; i++){
        (*b)(i) = (*b)(i) - a(i)/(*b)(i-1)*c(i-1);
        (*g)(i) = (*g)(i) - a(i)/(*b)(i-1)*(*g)(i-1);
    }
}

// Backwards substitution, general algorithm
void backward_gen(vec b, vec c, vec g, vec *v, int n){
    (*v)(n-1) = g(n-1)/b(n-1);
    for (int i = n-2; i >= 0; i--){
      (*v)(i) = (g(i) - c(i)*(*v)(i+1))/b(i);
    }
}

// Gaussian elimination, general algorithm
void gauss_elim_gen(vec a, vec b, vec c, vec g, vec *v, int n){
    forward_gen(a, &b, c, &g, n);
    backward_gen(b, c, g, v, n);
}



// Forward substitution, special algorithm
void forward_spec(vec a, vec *b, vec c, vec *g, int n){
    for (int i = 1; i < n; i++){
        (*b)(i) = 1.*(i+1)/i;
        (*g)(i) = (*g)(i) + (*g)(i-1)*(i-1)/i;
    }
}

// Backwards substitution, special algorithm
void backward_spec(vec b, vec c, vec g, vec *v, int n){
    (*v)(n-1) = g(n-1)/b(n-1);
    for (int i = n-2; i >= 0; i--){
      (*v)(i) = (g(i) + (*v)(i+1))/(i+1);
    }
}

// Gaussian elimination, special algorithm
void gauss_elim_spec(vec a, vec b, vec c, vec g, vec *v, int n){
    forward_spec(a, &b, c, &g, n);
    backward_spec(b, c, g, v, n);
}



// Writing to file exact_n.txt
void writetofile_exact(vec x, vec u, int n){
    ofstream myfile;
    string filename = "exact_" + to_string(n) + ".txt";
    myfile.open (filename);
    myfile << setw(15) << scientific << "x";
    myfile << setw(15) << scientific << "u" << endl;
    for (int i = 0; i < n; i++){
        myfile << setw(15) << scientific << x(i);
        myfile << setw(15) << scientific << u(i) << endl;
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
        myfile << setw(15) << scientific << v(i) << endl;
    }
    myfile.close();
}

// Writing to files exact_n.txt and approx_n.txt
void writetofile(vec x, vec u, vec v, int n){
    writetofile_exact(x, u, n);
    writetofile_approx(x, v, n);
}
