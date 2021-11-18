#include "functions.hpp"

using namespace arma;
using namespace std;


/**
 * inline function that implements periodic boundary conditions
 * returns the neighbouring indices of an element in the lattice
 *
 * i:       index of the element in question
 * limit:   length of the rows/columns of the square matrix
 * add:     number added to i
 */
inline int PBC(int i, int limit, int add){
    return (i + limit + add) % (limit);
}


/**
 * function that returns a square matrix of size L x L
 * L:       length of the rows/columns of the square matrix
 */
arma::mat spin_matrix(int L){
    mat S = mat(L, L);
    return S;
}


/**
 * function that initializes spin configuration, energy and magnetization
 * L:       length of the rows/columns of the square matrix
 * S:        of the rows/columns of the square matrix
 */

void initialize(int L, mat &S, double &E, double &M, int N, bool random){
    if (random){
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<int> distribution(0,1);
        for(int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){
                S(i, j) = distribution(gen);
                if (S(i,j) == 0){
                    S(i,j) += -1;       // random spin configuration with spins up and down
                }
                M += S(i, j);           // magnetization

            }
        }
    }

    else{
        for(int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){
                S(i,j) = 1.0;           // ordered spin configuration with spin up
                M += S(i, j);           // magnetization
            }
        }
    }


    for(int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
                E -= S(i, j) * S(PBC(i, L, -1), j) + S(i, j) * S(i, PBC(j, L, -1)); // energy


        }
    }

}


// function that returns the energy shift due to flipping a single spin
int delta_E(mat &S, int L, int i, int j){
    return 2*S(i,j)*(S(i, PBC(j, L, -1))
            + S(PBC(i, L, -1), j)
            + S(i, PBC(j, L, 1))
            + S(PBC(i, L, 1), j));
}

// metropolis algorithm
void metropolis(mat &S, int L, double T, double &E, double &M, int N_cycles, int N, bool write, int burnin){
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    vec boltzmann = zeros<mat>(17);

    // possible energies
    for(int de =-8; de <= 8; de+=4) boltzmann(de+8) = exp(-de/T);


    double E_exp = 0; double E_exp_sq = 0; double M_exp = 0; double M_exp_sq = 0;
    double eps_exp = 0; double eps_exp_sq = 0; double m_exp = 0; double m_exp_sq = 0;
    double cV = 0; double chi = 0; double e_exp = 0; double mag_exp = 0;

    vec epsilon_exp = vec(N_cycles);
    vec epsilon_samples = vec(N_cycles);
    vec magn_exp = vec(N_cycles);

    for (int i = 1; i <= burnin; i++){
        for(int j = 0; j < N; j++) {
            int x = distribution(gen)*L;
            int y = distribution(gen)*L;

            int dE = delta_E(S, L, x, y);

            if (dE <= 0){
                S(x, y) *= (-1);        // flips spin
            }

            else if (distribution(gen) <= boltzmann(dE+8) ){
                S(x,y) *= (-1);         // flips spin
            }

        }
    }


    for (int i = 0; i <= N_cycles; i++){
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


        /*
        // writing to file for problem 6
        epsilon_samples = E/N;
        ofstream eps_exp_file;
        eps_exp_file.open("eps_histo_1.bin", ios::binary | ios::app);
        eps_exp_file << setw(25) << epsilon_samples << endl;
        eps_exp_file.close();
        */



        E_exp += E;
        E_exp_sq += E*E;
        M_exp += abs(M);
        M_exp_sq += M*M;

        /*
        ofstream eps_exp_file;
        eps_exp_file.open("test.bin", ios::binary | ios::app);
        eps_exp_file << setw(25) << E_exp << endl;
        eps_exp_file.close();
        */

        double norm = 1./(((double) i)*N);
        epsilon_exp = E_exp*norm;
        magn_exp = M_exp*norm;

        /*

        ofstream eps_exp_file;
        eps_exp_file.open("m_exp_1_unordered.bin", ios::binary | ios::app);
        eps_exp_file << setw(25) << magn_exp << endl;
        eps_exp_file.close();
        */


        // for problem 5
        // double norm = 1./(((double) i)*N);
        // write this to file
        // epsilon_exp = E_exp*norm;

    }

    // problem 4
    e_exp = E_exp/N;
    eps_exp = e_exp/N_cycles;
    //eps_exp_sq = E_exp_sq/(N * N * N_cycles);

    mag_exp = M_exp/N;
    m_exp = mag_exp/N_cycles;
    //m_exp_sq = M_exp_sq/(N * N * N_cycles);

    cV = 1./(T*T)*(E_exp_sq/N_cycles - E_exp/N_cycles * E_exp/N_cycles)/N;
    chi = 1./T*(M_exp_sq/N_cycles - M_exp/N_cycles * M_exp/N_cycles)/N;


    /*
    cout << setprecision(15) << eps_exp << endl;
    cout << setprecision(15) << m_exp << endl;
    cout << setprecision(15) << cV << endl;
    cout << setprecision(15) << chi << endl;
    */



    if (write){
        ofstream eps_exp_file;
        eps_exp_file.open("eps_exp_" + to_string(L) + "L.bin", ios::binary | ios::app);
        eps_exp_file << setw(25) << eps_exp << " " << T << endl;
        eps_exp_file.close();

        ofstream m_exp_file;
        m_exp_file.open("m_exp_" + to_string(L) + "L.bin", ios::binary | ios::app);
        m_exp_file << setw(25) << m_exp << " " << T << endl;
        m_exp_file.close();

        ofstream cV_file;
        cV_file.open("cV_" + to_string(L) + "L.bin", ios::binary | ios::app);
        cV_file << setw(25) << cV << " " << T << endl;
        cV_file.close();

        ofstream chi_file;
        chi_file.open("chi_" + to_string(L) + "L.bin", ios::binary | ios::app);
        chi_file << setw(25) << chi << " " << T << endl;
        chi_file.close();
    }
}
