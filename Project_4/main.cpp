#include "functions.hpp"

using namespace arma;
using namespace std;


int main(int argc, char const *argv[]){
    int L = atoi(argv[1]);               // lattice length
    double N = L*L;                      // number of spins
    double T = 2.1;
    int NT = 100;
    vec T_vec = linspace(2.1, 2.4, NT);  // cV peak between 2.25 and 2.35, the others peak at 2.1??
    int N_cycles = pow(2, 21);


    //mat S = spin_matrix(L);

    double E = 0.;                       // initialize energy
    double M = 0.;                       // initialize magnetization

    //initialize(L, S, E, M, N);
    //metropolis(S, L, T, E, M, N_cycles, N);



    bool timing = true;
    if (timing){
        auto t0 = omp_get_wtime();
        #pragma omp parallel private(E, M) // Start parallel region
        {
            mat S = spin_matrix(L);
            initialize(L, S, E, M, N);
            #pragma omp for
            for(int i = 0; i < NT; i++){
                metropolis(S, L, T_vec(i), E, M, N_cycles, N);
            }
        }
        auto t1 = omp_get_wtime();
        double duration_seconds_wo = std::chrono::duration<double>(t1 - t0).count();

        cout << "time used with OpenMP = " << duration_seconds_wo << " seconds\n";
    }



    return 0;
}
