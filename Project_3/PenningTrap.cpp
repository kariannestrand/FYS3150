#include "PenningTrap.hpp"
#include "Particle.hpp"

using namespace arma;
using namespace std;


PenningTrap::PenningTrap(double B0, double V0, double d, double ke, double f, vec omega_v, int n, double N, mat pos, mat vel, vec q_vec, vec m_vec, bool write, bool interaction, bool modified){
    B0_ = B0;                   // magnetic field strength
    V0_ = V0;                   // applied potential
    d_ = d;                     // characteristic dimension
    ke_ = ke;
    f_ = f;
    omega_v_ = omega_v;
    n_ = n;
    N_ = N;
    write_ = write;
    interaction_ = interaction;
    modified_ = modified;

    // making list/contatiner for particle objects

    for (int i = 0; i < n_; i++){
        particles_.push_back(Particle(q_vec(i), m_vec(i), pos.col(i), vel.col(i)));
    }

}


vec PenningTrap::external_B_field(int i){
    Particle& p_i = particles_[i];
    vec r = p_i.r_;

    vec B = vec(3).fill(0.);

    if (modified_){
        // constraints on magnetic field
        if (norm(r) <= d_){
            B(2) = B0_;
        }
    }
    else{
        B(2) = B0_;
    }


    return B;
}


vec PenningTrap::external_E_field(int i, int k, double dt){
    Particle& p_i = particles_[i];
    vec r = p_i.r_;

    vec F = vec(3).fill(0);
    // position of E-field
    F(0) = -1.;
    F(1) = -1.;
    F(2) = 2.;

    vec E = vec(3).fill(0);

    if (modified_){
        // constraints on electric field
        if (norm(r) <= d_){
            double V0t = V0_*(1 + f_*cos(omega_v_(k)*i*dt));
            E = - V0t/(d_*d_)*F % r;
        }
    }
    else{
        E = - V0_/(d_*d_)*F % r;
    }

    return E;
}


vec PenningTrap::force_particle(int i, int j){
    Particle& p_i = particles_[i];
    Particle& p_j = particles_[j];
    vec r = p_i.r_ - p_j.r_;

    int q_i = p_i.q_;
    int q_j = p_j.q_;

    double dr = norm(r) * norm(r) * norm(r);


    vec F = vec(3).fill(0);
    F = ke_*(q_i*q_j)/dr * r;

    return F;
}


vec PenningTrap::total_force_external(int i, int k, double dt){
    Particle& p_i = particles_[i];

    vec F = vec(3).fill(0);
    vec E = external_E_field(i, k, dt);
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


vec PenningTrap::total_force(int i, int k, double dt){
    vec F = vec(3).fill(0);

    if (interaction_){
        F = total_force_external(i, k, dt) + total_force_particles(i);
    }
    else{
        F = total_force_external(i, k, dt);
    }

    return F;

}


int PenningTrap::particles_trapped(){
    int count = 0;
    for (int i=0; i < particles_.size(); i++){
        vec r = particles_[i].r_;
        if (norm(r) <= d_){
        count = count + 1;
        }
    }

    return count;
}


void PenningTrap::evolve_RK4(double dt, int k){
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
            vec F = total_force(i, k, dt);
            vec a = F/p_i.m_;

            vec K1_v = a*dt;
            vec K1_r = p_i.v_*dt;


            // K2
            p_i.v_ = v_old + (1/2.)*K1_v;
            p_i.r_ = r_old + (1/2.)*K1_r;

            F = total_force(i, k, 0.5*dt);
            a = F/p_i.m_;

            vec K2_v = a*dt;
            vec K2_r = p_i.v_*dt;


            // K3
            p_i.v_ = v_old + (1/2.)*K2_v;
            p_i.r_ = r_old + (1/2.)*K2_r;

            F = total_force(i, k, 0.5*dt);
            a = F/p_i.m_;

            vec K3_v = a*dt;
            vec K3_r = p_i.v_*dt;


            // K4
            p_i.v_ = v_old + (1/2.)*K3_v;
            p_i.r_ = r_old + (1/2.)*K3_r;

            F = total_force(i, k, dt);
            a = F/p_i.m_;

            vec K4_v = a*dt;
            vec K4_r = p_i.v_*dt;

            // last step
            p_i.v_ = v_old + (1/6.)*(K1_v + 2.*K2_v + 2.*K3_v + K4_v);
            p_i.r_ = r_old + (1/6.)*(K1_r + 2.*K2_r + 2.*K3_r + K4_r);

            V_total.slice(j).col(i) = p_i.v_;
            R_total.slice(j).col(i) = p_i.r_;

        }
    }

    if (modified_){
        int N = particles_trapped();

        ofstream file;
        file.open("trapped_0001dt_002w_01f.txt", ios::app);
        file << setw(25) << N << " " << omega_v_(k) << endl;
        file.close();
    }


    if (write_){
        for (int i = 0; i < n_; i++){
            mat R = mat(3, n_);
            mat V = mat(3, n_);
            R = R_total.col(i);
            V = V_total.col(i);
            R.save("RK4_r_" + to_string(i) + "_0001dt" + ".bin");
            V.save("RK4_v_" + to_string(i) + "_0001dt" + ".bin");
        }
    }

}



void PenningTrap::evolve_forward_Euler(double dt, int k){
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

            vec F = total_force(i, k, dt);
            vec a = F/p_i.m_;


            p_i.v_ = p_i.v_ + a*dt;
            p_i.r_ = p_i.r_ + p_i.v_*dt;


            V_total.slice(j).col(i) = p_i.v_;
            R_total.slice(j).col(i) = p_i.r_;
        }
    }


    if (write_){
        for (int i = 0; i < n_; i++){
            mat R = mat(3, n_);
            mat V = mat(3, n_);
            R = R_total.col(i);
            V = V_total.col(i);
            R.save("Euler_r_" + to_string(i) + "_0001dt" + ".bin");
            V.save("Euler_v_" + to_string(i) + "_0001dt" + ".bin");
        }
    }

}
