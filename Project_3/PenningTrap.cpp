
#include "PenningTrap.hpp"
#include "Particle.hpp"

using namespace arma;
using namespace std;


PenningTrap::PenningTrap(double B0, double V0, double d, double ke, int n, double N, mat pos, mat vel, vec q_vec, vec m_vec, bool write, bool interaction){
    B0_ = B0;                   // magnetic field strength
    V0_ = V0;                   // applied potential
    d_ = d;                     // characteristic dimension
    ke_ = ke;
    n_ = n;
    N_ = N;
    write_ = write;
    interaction_ = interaction;

    // making list/contatiner for particle objects

    for (int i = 0; i < n_; i++){
        particles_.push_back(Particle(q_vec(i), m_vec(i), pos.col(i), vel.col(i)));
    }

}


vec PenningTrap::external_B_field(int i){
    Particle& p_i = particles_[i];
    vec r = p_i.r_;

    vec B = vec(3).fill(0.);

    B(2) = B0_;

    /*
    if ((r(0) > d_) && (r(1) > d_) && (r(2) > d_)){
        B(2) = B0_;
    }
    */

    return B;
}


vec PenningTrap::external_E_field(int i){
    Particle& p_i = particles_[i];
    vec r = p_i.r_;

    vec F = vec(3).fill(0);
    // position of E-field
    F(0) = -1.;
    F(1) = -1.;
    F(2) = 2.;

    vec E = - V0_/(d_*d_)*F % r;

    /*
    vec E;
    if ((r(0) > d_) && (r(1) > d_) && (r(2) > d_)){
        E = - V0_/(d_*d_)*F % r;
    }
    else{
        E(0) = 0;
        E(1) = 0;
        E(2) = 0;
    }
    */
    return E;
}


vec PenningTrap::force_particle(int i, int j){
    Particle& p_i = particles_[i];
    Particle& p_j = particles_[j];
    vec r = p_i.r_ - p_j.r_;

    int q_i = p_i.q_;
    int q_j = p_j.q_;

    vec dr = abs(r) % abs(r) % abs(r);


    vec F;
    if (interaction_){
        F = ke_*(q_i*q_j)/dr % r;
    }
    else{
        F(0) = 0;
        F(1) = 0;
        F(2) = 0;
    }


    return F;
}


vec PenningTrap::total_force_external(int i){
    Particle& p_i = particles_[i];

    vec F = vec(3).fill(0);
    vec E = external_E_field(i);
    vec B = external_B_field(i);

    vec v = p_i.v_;
    int q = p_i.q_;

    F = q*E + cross(q*v, B);

    return F;
}


vec PenningTrap::total_force_particles(int i){
    vec F = vec(3).fill(0);

    for (int j = 0; j < n_; j++){
        if (j != i){
            F = F + force_particle(i, j);
        }
    }

    return F;
}


vec PenningTrap::total_force(int i){
    vec F = vec(3).fill(0);

    if (interaction_){
        F = total_force_external(i) + total_force_particles(i);
    }
    else{
        F = total_force_external(i);
    }

    return F;

}


void PenningTrap::evolve_RK4(double dt){
    cube R_total = cube(3, n_, N_);
    cube V_total = cube(3, n_, N_);


    for (int i = 0; i < n_; i++){
        Particle& p_i = particles_[i];
        R_total.slice(0).col(i) = p_i.r_;
        V_total.slice(0).col(i) = p_i.v_;
    }


    for (int j = 1; j < N_; j++){

        for (int i = 0; i < n_; i++){
            Particle& p_i = particles_[i];

            vec r_old = R_total.slice(j-1).col(i);
            vec v_old = V_total.slice(j-1).col(i);

            // K1
            vec F = total_force(i);
            vec a = F/p_i.m_;

            vec K1_v = a*dt;
            vec K1_r = p_i.v_*dt;


            // K2
            p_i.v_ = v_old + (1/2.)*K1_v;
            p_i.r_ = r_old + (1/2.)*K1_r;

            F = total_force(i);
            a = F/p_i.m_;

            vec K2_v = a*dt;
            vec K2_r = p_i.v_*dt;


            // K3
            p_i.v_ = v_old + (1/2.)*K2_v;
            p_i.r_ = r_old + (1/2.)*K2_r;

            F = total_force(i);
            a = F/p_i.m_;

            vec K3_v = a*dt;
            vec K3_r = p_i.v_*dt;


            // K4
            p_i.v_ = v_old + (1/2.)*K3_v;
            p_i.r_ = r_old + (1/2.)*K3_r;

            F = total_force(i);
            a = F/p_i.m_;

            vec K4_v = a*dt;
            vec K4_r = p_i.v_*dt;

            // last step
            p_i.v_ = v_old + (1/6.)*(K1_v + 2.*K2_v + 2.*K3_v + K4_v);
            p_i.r_ = r_old + (1/6.)*(K1_r + 2.*K2_r + 2.*K3_r + K4_r);

            V_total.slice(j).col(i) = p_i.v_;
            R_total.slice(j).col(i) = p_i.r_;

            // p_i.r_ = r_old;
            // p_i.v_ = v_old;

            /*
            if (i == 0){
                if (write_){
                    ofstream file;
                    file.open("RK4_v_1_0001dt.txt", ios::app);
                    file << V(0, i) << " " << V(1, i) << " " << V(2, i) << endl;
                    file.close();
                }

                if (write_){
                    ofstream file;
                    file.open("RK4_r_1_0001dt.txt", ios::app);
                    file << R(0, i) << " " << R(1, i) << " " << R(2, i) << endl;
                    file.close();
                }
            }
            if (i == 1){
                if (write_){
                    ofstream file;
                    file.open("RK4_v_2_0001dt.txt", ios::app);
                    file << V(0, i) << " " << V(1, i) << " " << V(2, i) << endl;
                    file.close();
                }

                if (write_){
                    ofstream file;
                    file.open("RK4_r_2_0001dt.txt", ios::app);
                    file << R(0, i) << " " << R(1, i) << " " << R(2, i) << endl;
                    file.close();
                }
            }
            */

        }


    }


    if (write_){
        for (int i = 0; i < n_; i++){
            mat R = mat(3, n_);
            mat V = mat(3, n_);
            R = R_total.col(i);
            V = V_total.col(i);
            R.save("r_" + to_string(i) + "_0001" + ".bin");
            V.save("v_" + to_string(i) + "_0001" + ".bin");

        }
    }
}


void PenningTrap::evolve_forward_Euler(double dt){
    mat R = mat(3, n_).fill(0);
    mat V = mat(3, n_).fill(0);

    for (int i = 0; i < n_; i++){
        Particle& p_i = particles_[i];


        for (int j = 0; j < N_; j++){

            vec r_old = p_i.r_;
            vec v_old = p_i.v_;

            vec F = total_force(i);
            vec a = F/p_i.m_;


            p_i.v_ = p_i.v_ + a*dt;
            p_i.r_ = p_i.r_ + p_i.v_*dt;


            V.col(i) = p_i.v_;
            R.col(i) = p_i.r_;

            if (i == 0){
                if (write_){
                    ofstream file;
                    file.open("Euler_v_1_0001dt.txt", ios::app);
                    file << V(0, i) << " " << V(1, i) << " " << V(2, i) << endl;
                    file.close();
                }

                if (write_){
                    ofstream file;
                    file.open("Euler_r_1_0001dt.txt", ios::app);
                    file << R(0, i) << " " << R(1, i) << " " << R(2, i) << endl;
                    file.close();
                }
            }
            else if (i == 1){
                if (write_){
                    ofstream file;
                    file.open("Euler_v_2_0001dt.txt", ios::app);
                    file << V(0, i) << " " << V(1, i) << " " << V(2, i) << endl;
                    file.close();
                }

                if (write_){
                    ofstream file;
                    file.open("Euler_r_2_0001dt.txt", ios::app);
                    file << R(0, i) << " " << R(1, i) << " " << R(2, i) << endl;
                    file.close();
                }
            }
        }
    }
}
