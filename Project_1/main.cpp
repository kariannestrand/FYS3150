#include "functions.hpp"

using namespace arma;
using namespace std;


int main(int argc, char* argv[]){
    int n = atoi(argv[1]);                 // length of array
    double h = 1./n;                       // step size
    vec a = vec(n).fill(-1.);              // defining a-vector filled with -1's
    vec b_gen = vec(n).fill(2.);           // defining b-vector filled with 2's for general algortihm
    vec b_spec = vec(n).fill(2.);          // defining b-vector filled with 2's for special algortihm
    vec c = vec(n).fill(-1.);              // defining c-vector filled with -1's
    vec x = arma::linspace(0+h, 1-h, n);   // defining array x in [0, 1]
    vec f = 100*exp(-10*x);                // defining source term
    vec g_gen = h*h*f;                     // defining solution vector g for general algortihm
    vec g_spec = h*h*f;                    // defining solution vector g for special algortihm
    vec v_gen = vec(n).fill(0.);           // creating empty solution vector v for general algortihm
    vec v_spec = vec(n).fill(0.);          // creating empty solution vector v gor special algorithm

    vec u = exact(x);                      // calling exact solution function


    // Measures duration of general and special algorithm and prints outcome if true
    bool timing = true;
    if (timing){
        auto t1_gen = std::chrono::high_resolution_clock::now();
        gauss_elim_gen(a, b_gen, c, g_gen, &v_gen, n);
        auto t2_gen = std::chrono::high_resolution_clock::now();
        double duration_seconds_gen = std::chrono::duration<double>(t2_gen - t1_gen).count();

        auto t1_spec = std::chrono::high_resolution_clock::now();
        gauss_elim_spec(a, b_spec, c, g_spec, &v_spec, n);
        auto t2_spec = std::chrono::high_resolution_clock::now();
        double duration_seconds_spec = std::chrono::duration<double>(t2_spec - t1_spec).count();

        cout << "time used by the general algorithm = " << duration_seconds_gen << " seconds\n";
        cout << "time used by the special algorithm = " << duration_seconds_spec << " seconds\n";
    }
    else{
            gauss_elim_gen(a, b_gen, c, g_gen, &v_gen, n);
            gauss_elim_spec(a, b_spec, c, g_spec, &v_spec, n);
    }

    // Creates files exact_n.txt and approx_n.txt if true
    bool write_to_file = false;
    if (write_to_file){
        writetofile(x, u, v_gen, n);
    }

    return 0;
}
