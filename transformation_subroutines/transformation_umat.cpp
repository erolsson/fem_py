//
// Created by erolsson on 11/09/2019.
//

#include "transformation_umat.h"

#include <fstream>
#include <iostream>
#include <vector>

#include "Eigen/Dense"

#include "simulation_parameters.h"
#include "stress_functions.h"
#include "transformation_umat.h"

double yield_function(const Eigen::Matrix<double, 6, 1>& stilde, double sigma_y) {
    return 3*double_contract(stilde, stilde)/2 - sigma_y*sigma_y;
}
/*
double transformation_function(const Eigen::Matrix<double, 6, 1>& stress, double fm, double k, double Ms, double Mss,
        double T, std::array<double, 3> a) {
    return 1 - exp(k*(Ms + ms_stress(stress, )))
}

double ms_stress(const Eigen::Matrix<double, 6, 1>& stress, double a1, double a2, double a3) {
    Eigen::Matrix<double, 6, 1> s_dev = deviator(stress);
    double m_stress = a1*(stress[0] + stress[1] + stress[2]);       // Contribution from hydrostic stress
    m_stress += a2*von_Mises(stress);
    m_stress += a3*vector_det(s_dev);
    return m_stress;
}
*/

class State {
public:
    using Vector6 = Eigen::Matrix<double, 6, 1>;
    explicit State(double* data, unsigned back_stresses) : data_(data), back_stresses_(back_stresses) {}
    double& ep_eff() { return data_ [0]; }
    double& fM() { return  data_[1]; }
    double& R() { return data_ [2]; }

    Eigen::Map<Vector6> back_stress_vector(unsigned n) {
        return Eigen::Map<Vector6>(&data_[3 + n*6]);
    }

    Eigen::Map<Vector6> total_back_stress() {
        return Eigen::Map<Vector6>(&data_[3 + back_stresses_*6]);
    }

private:
    double* data_;
    unsigned back_stresses_;
};

extern "C" void umat_(double *stress, double *statev, double *ddsdde, double *sse, double *spd, double *scd,
        double *rpl, double *ddsddt, double *drplde, double *drpldt, double *stran, double *dstran, double *time,
        double *dtime, double *temp, double *dtemp, double *predef, double *dpred, char *cmname, int *ndi, int *nshr,
        int *ntens, int *nstatv, const double* props, int *nprops, double *coords, double *drot, double *pnewdt,
        double *celent, double *dfgrd0, double *dfgrd1, int *noel, int *npt, int *layer, int *kspt, int *kstep,
        int *kinc, short cmname_len) {

    using Matrix6x6 = Eigen::Matrix<double, 6, 6, Eigen::RowMajor>;
    using Vector6 = Eigen::Matrix<double, 6, 1>;

    TransformationMaterialParameters params(props);
    double G = params.E()/2/(1+params.v());
    double K = params.E()/3/(1-2*params.v());
    // Collecting state variables
    State state(statev, params.back_stresses());

    double sy = params.sy0() + state.R();
    // Vector6 back_stress;

    Eigen::Map<Matrix6x6> D_alg(ddsdde);
    Matrix6x6 Del = 2*G*J + K*E3;

    Eigen::Map<Vector6> stress_vec(stress);
    Eigen::Map<Vector6> de(dstran);

    Vector6  st = stress_vec + Del*de;  // Trial stress

    Vector6 stilde = deviator(st);
    if (params.kinematic_hardening()) {
        stilde -= state.total_back_stress();
    }
    bool plastic = params.plastic() && yield_function(stilde, sy) > 0;
    bool phase_transformations = params.strain_transformation || params.stress_transformation;
    phase_transformations = false;
    bool elastic = !plastic && !phase_transformations;
    stress_vec = st;
    std::cout << "elastic:" << elastic << std::endl;
    if (elastic) {     // Use the trial stress as the stress and the elastic stiffness matrix as the tangent
        D_alg = Del;
    }
    else {
        double DfM = 0;
        double RA = 0;

        Vector6 nij;
        Vector6 st_dev = deviator(st);

        Matrix6x6 A = I;   // Used to formulate the tangent matrix DDSDDE = A^-1 * Del


        // Calculating the increment in effective plastic strain DL
        double DL = 0;
        double residual = 1e99;

        double ds2eqdDL = -3*G;
        double dsydDL = 0;

        double D = 0;
        double R2 = 0;
        Eigen::VectorXd theta = Eigen::VectorXd::Zero(params.back_stresses());
        Eigen::VectorXd Am = Eigen::VectorXd::Zero(params.back_stresses());
        while(residual > 1e-15) {
            Vector6 s_prime = st_dev;

            if (params.isotropic_hardening()) {
                R2 = (state.R() + params.b()*params.Q()*DL)/(1 + params.b()*DL);
                sy = params.sy0() + R2;
                D += R2;
                dsydDL = params.b()*(params.Q() - state.R());
            }

            double Cm_theta_m = 0;
            if (params.kinematic_hardening()) {
                for (unsigned i=0; i != params.back_stresses(); ++i) {
                    double Cm = params.Cm(i);
                    double gm = params.gamma(i);

                    theta[i] = 1./(1+gm*DL);
                    Am[i] = gm/((1+gm*DL)*(1+gm*DL));

                    s_prime -= theta[i]*state.back_stress_vector(i);
                    Cm_theta_m += theta[i]*Cm;
                }
            }

            double s_eq_prime = sqrt(1.5*double_contract(s_prime, s_prime));
            nij = 1.5*s_prime/s_eq_prime;

            if (params.kinematic_hardening()) {
                Vector6 ds_prime_dDL = Vector6::Zero();
                for (unsigned i=0; i != params.back_stresses(); ++i) {
                    ds_prime_dDL += Am[i]*state.back_stress_vector(i);
                    ds2eqdDL -= theta[i]*theta[i]*params.Cm(i)*DL;
                }
                ds2eqdDL += double_contract(nij, ds_prime_dDL);
            }
            double s_eq_2 = s_eq_prime - 3*G*(DL+R2*DfM) - Cm_theta_m*DL;
            if (phase_transformations) {
                double numerator = 1+3*G*R2*DfM/params.sy0_au();
                s_eq_2 /= numerator;
                ds2eqdDL /= numerator;
            }

            double f = s_eq_2 - sy;

            double dfdDL = ds2eqdDL - dsydDL;
            if (! phase_transformations) {
                double dDL = f/dfdDL;
                DL -= dDL;
                residual = abs(dDL);
            }

            std::cout << "DL=" << DL << "  s_eq_2=" << s_eq_2 << " s_eq_prime=" << s_eq_prime << std::endl;
        }

        state.ep_eff() += DL;
        state.R() = R2;
        stress_vec -= 2*G*DL*nij;

        // Helpful quantity in the calculation of the tangent DDSDDE
        double Y = sy;
        double H = 0;
        double Csum = 0;
        double thetaCsum = 0;
        if (params.isotropic_hardening()) {
            H += dsydDL;
        }
        // Updating back stresses
        if (params.kinematic_hardening()) {
            state.total_back_stress() *= 0;
            Vector6 sum_back_stress_gamma = Vector6::Zero();
            for (unsigned i = 0; i != params.back_stresses(); ++i) {
                state.back_stress_vector(i) += 2./3*params.Cm(i)*DL*nij;
                state.back_stress_vector(i) *= theta[i];
                state.total_back_stress() += state.back_stress_vector(i);
                Y += theta[i]*params.Cm(i)*DL;
                sum_back_stress_gamma += state.back_stress_vector(i)*params.gamma(i);
                Csum += params.Cm(i);
                thetaCsum += params.Cm(i)*theta(i);
            }
            H += Csum - double_contract(nij, sum_back_stress_gamma);
        }

        // Calculating the tangent modulus
        // Ideal plasticity, i. e no hardening is a special case as the derivation assumes H != 0
        Matrix6x6 nnt = nij*nij.transpose();
        if (H != 0) {
            A += 3*G/Y*(DL+RA*DfM)*J;
            A += 2*G*(1./H - (DL + RA*DfM)/Y/H*(thetaCsum + dsydDL))*nnt;
            if (params.kinematic_hardening()) {
                Vector6 adjusted_back_stress = Vector6::Zero();
                for (unsigned i = 0; i != params.back_stresses(); ++i) {
                    adjusted_back_stress += Am[i]/theta[i]*state.back_stress_vector(i);
                }
                A += 3*G/H/Y*(DL + RA*DfM)*adjusted_back_stress*nij.transpose();
            }
        }
        A *= 2;
        A.block(0, 0, 3, 3) /=2;
        D_alg = A.inverse()*Del;
        // D_alg = Del;

    }
}