#include <armadillo>
#include <iostream>

using namespace arma;
using namespace std;


arma::mat spin_matrix(int L);
void initialize(int L, mat &S, double &E, double &M, int N);
int delta_E(mat &S, int L, int i, int j);
void metropolis(mat &S, int L, double T, double &E, double &M, int N_cycles, int N);

// inline function for Periodic Boundary Conditions
inline int PBC(int i, int limit, int add){
    return (i + limit + add) % (limit);
}


int main(int argc, char const *argv[]){
    int L = 20;                         // lattice length
    double N = L*L;                     // number of spins
    double T = 1.0;
    int N_cycles = 1000;
    

    mat S = spin_matrix(L);

    double E = 0.;                      // initialize energy
    double M = 0.;                      // initialize magnetization

    initialize(L, S, E, M, N);
    metropolis(S, L, T, E, M, N_cycles, N);


    return 0;
}


// matrix for number of spins
arma::mat spin_matrix(int L){
    mat S = mat(L,L);
    return S;
}


// function to initialize spin configuration, energy and magnetization
void initialize(int L, mat &S, double &E, double &M, int N){
    bool random = false;        // random if true, ordered if false

    if (random){
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<int> distribution(0,1);
        for(int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){
                S(i, j) = distribution(gen);        
                
                if (S(i,j) == 0){
                    S(i,j) += -1;
                }
                M += S(i, j);

            }
        }
    }
    
    else{
        for(int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){
                S(i,j) = 1.0;           // spin configuration, choose -1 for down, 1 for up
                M += S(i, j);
            }
        }
    }


    for(int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
                E -= S(i, j) * S(PBC(i, L, -1), j) + S(i, j) * S(i, PBC(j, L, -1));
                

        }
    }
    
}

// delta E
int delta_E(mat &S, int L, int i, int j){
    return 2*S(i,j)*(S(i, PBC(j,L,-1))
            + S(PBC(i,L,-1),j) 
            + S(i, PBC(j,L,1)) 
            + S(PBC(i,L,1),j)); 
}

// metropolis algorithm
void metropolis(mat &S, int L, double T, double &E, double &M, int N_cycles, int N){

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    //made from mortens lecture notes
    vec boltzmann = zeros<mat>(17);
    
    // possible energies
    for(int de =-8; de <= 8; de+=4) boltzmann(de+8) = exp(-de/T);
    

    double E_exp = 0; double E_exp_sq = 0; double M_exp = 0; double M_exp_sq = 0;
    double eps_exp 0; double eps_exp_sq 0; double m_exp 0; double m_exp_sq 0;

    vec epsilon_exp = vec(N_cycles); 


    for (int i = 1; i <= N_cycles; i++){
        for(int j = 0; j < N; j++) {
            int x = distribution(gen)*L;
            int y = distribution(gen)*L;

            int dE = delta_E(S, L, x, y);

            if (dE <= 0){
                S(x, y) *= (-1);        // flips spin
                E += dE;
                M += 2*S(x, y);
            }
            
            else if (distribution(gen) <= boltzmann(dE+8) ){
                S(x,y) *= (-1);         // flips spin
                E += dE;
                M += 2*S(x, y);
                
            }

        }


        // for histogram in problem 6
        // write this to file
        epsilon_samples = E/(N*N);

        E_exp += E;
        E_exp_sq += E*E;
        M_exp += abs(M);
        M_exp_sq += M*M;
        
        
        // for problem 5
        double norm = 1./(((double) i)*N);
        // write this to file
        epsilon_exp = E_exp*norm;
        
        
        
    
    }


    // problem 4
    eps_exp /= N * N_cycles;
    eps_exp_sq /= N * N * N_cycles;
    m_exp /= N * N_cycles;
    m_exp_sq /= N * N * N_cycles;
    
    /*
    cout << eps_exp << endl;
    cout << eps_exp_sq << endl;
    cout << m_exp << endl;
    cout << m_exp_sq << endl;
    */
    
}


