#include "functions.hpp"

using namespace arma;
using namespace std;


int main(int argc, char const *argv[]){
    bool random = false;                 // random spin configuration if true, ordered spin configuration if false
    bool one_temperature = false;        // runs program for one temperature (second input argument after ./main.exe) if true, multiple temperatures (T_vec) if false
    bool write = true;                   // writes bin-files if true, does not if false

    int L = atoi(argv[1]);               // lattice length
    double N = L*L;                      // number of spins
    int NT = 100;                        // number of temperature steps
    vec T_vec = linspace(2.1, 2.4, NT);  // temperature vector
    int N_cycles = pow(2, 21);           // number of cycles
    int burnin = 50000;

    double E = 0.;                       // initialize energy
    double M = 0.;                       // initialize magnetization


    if (one_temperature){
        auto t0 = std::chrono::high_resolution_clock::now();   // start clock

        double T = atoi(argv[2]);
        mat S = spin_matrix(L);
        initialize(L, S, E, M, N, random);
        metropolis(S, L, T, E, M, N_cycles, N, write, burnin);

        auto t1 = std::chrono::high_resolution_clock::now();   // end clock
        double duration_seconds = std::chrono::duration<double>(t1 - t0).count();                                // calculates duration
        cout << "Time used with one temperature (no parallelization) = " << duration_seconds << " seconds\n";    // prints duration
    }
    else{
        auto t0 = omp_get_wtime();   // start clock

        #pragma omp parallel private(E, M)  // Start parallel region
        {
            mat S = spin_matrix(L);
            initialize(L, S, E, M, N, random);
            #pragma omp for
            for(int i = 0; i < NT; i++){
                metropolis(S, L, T_vec(i), E, M, N_cycles, N, write, burnin);
            }
        }

        auto t1 = omp_get_wtime();   // end clock
        double duration_seconds = std::chrono::duration<double>(t1 - t0).count();                   // calculates duration
        cout << "Time used with OpenMP parallelization = " << duration_seconds << " seconds\n";     // prints duration
    }

    return 0;
}
