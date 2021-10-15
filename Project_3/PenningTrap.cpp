#include "PenningTrap.hpp"
#include "Particle.hpp"

using namespace arma;
using namespace std;

PenningTrap::PenningTrap(double B0, double V0, double d, double ke, int n, double N, mat pos, mat vel, vec q_vec, vec m_vec){
    B0_ = B0;                   // magnetic field strength
    V0_ = V0;                   // applied potential
    d_ = d;                     // characteristic dimension
    ke_ = ke;
    n_ = n;
    N_ = N;

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
    //if ((r(0) > d_) && (r(1) > d_) && (r(2) > d_)){
    E = - V0_/(d_*d_)*F % r;
    //}
    /*
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

    vec F = ke_*(q_i*q_j)/dr % r;

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
    F = total_force_external(i) + total_force_particles(i);

    return F;

}


void PenningTrap::evolve_RK4(double dt, bool write){
    mat R = mat(3, n_).fill(0);
    mat V = mat(3, n_).fill(0);

    for (int j = 0; j < N_; j++){
        for (int i = 0; i < n_; i++){
            Particle& p_i = particles_[i];

            // K1
            vec F = total_force(i);
            vec a = F/p_i.m_;

            vec K1_v = a*dt;
            vec K1_r = p_i.v_*dt;

            p_i.v_ = p_i.v_ + (1/2.)*K1_v;
            p_i.r_ = p_i.r_ + (1/2.)*K1_r;


            // K2
            F = total_force(i);
            a = F/p_i.m_;

            vec K2_v = p_i.v_ + (1/2.)*K1_v;
            vec K2_r = p_i.r_ + (1/2.)*K1_r;

            p_i.v_ = p_i.v_ + (1/2.)*K2_v;
            p_i.r_ = p_i.r_ + (1/2.)*K2_r;


            // K3
            F = total_force(i);
            a = F/p_i.m_;

            vec K3_v = p_i.v_ + (1/2.)*K2_v;
            vec K3_r = p_i.r_ + (1/2.)*K2_r;

            p_i.v_ = p_i.v_ + (1/2.)*K3_v;
            p_i.r_ = p_i.r_ + (1/2.)*K3_r;

            // K4
            F = total_force(i);
            a = F/p_i.m_;

            vec K4_v = p_i.v_ + K3_v;
            vec K4_r = p_i.r_ + K3_r;

            p_i.v_ = p_i.v_ + K4_v;
            p_i.r_ = p_i.r_ + K4_r;

            // last step
            V.col(i) = p_i.v_ + (1/6.)*(K1_v + 2.*K2_v + 2.*K3_v + K4_v);
            R.col(i) = p_i.r_ + (1/6.)*(K1_r + 2.*K2_r + 2.*K3_r + K4_r);

            

        }

        if (write){
            ofstream file;
            file.open("RK4_v.txt", ios::app);
            file << V << endl;
            file.close();
        }

        if (write){
            ofstream file;
            file.open("RK4_r.txt", ios::app);
            file << R << endl;
            file.close();
        }
    }

}


void PenningTrap::evolve_forward_Euler(double dt, bool write){
    mat R = mat(3, n_).fill(0);
    mat V = mat(3, n_).fill(0);
    for (int j = 0; j < N_; j++){
        for (int i = 0; i < n_; i++){
            Particle& p_i = particles_[i];

            //cout << p_i.r_ << endl;
            vec F = total_force(i);
            vec a = F/p_i.m_;


            p_i.v_ = p_i.v_ + a*dt;
            p_i.r_ = p_i.r_ + p_i.v_*dt;


            V.col(i) = p_i.v_;
            R.col(i) = p_i.r_;
        }

        if (write){
            ofstream file;
            file.open("Euler_v.txt", ios::app);
            file << V << endl;
            file.close();
        }

        if (write){
            ofstream file;
            file.open("Euler_r.txt", ios::app);
            file << R << endl;
            file.close();
        }
    }

}
