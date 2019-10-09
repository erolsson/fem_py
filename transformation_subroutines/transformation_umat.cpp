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
        return Eigen::Map<Vector6>(data_ + 3 + (n-1)*6);
    }

    Eigen::Map<Vector6> total_back_stress(unsigned n) {
        return Eigen::Map<Vector6>(data_ + 3 + back_stresses_*6);
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
    bool plastic = params.plastic() && yield_function(stilde, sy) > 0;
    bool phase_transformations = params.strain_transformation || params.stress_transformation;
    phase_transformations = false;
    bool elastic = !plastic && !phase_transformations;
    stress_vec = st;
    if (elastic) {     // Use the trial stress as the stress and the elastic stiffness matrix as the tangent
        D_alg = Del;
    }
    else {
        double DfM = 0;
        double RA = 0;

        Vector6 nij;
        Vector6 st_dev = deviator(st);
        if (plastic) {
            // Calculating the increment in effective plastic strain DL
            double DL = 0;
            double dDL = 1e99;

            double dsydDL = 0;

            double D = 0;
            double R2 = 0;
            Eigen::VectorXd theta = Eigen::VectorXd::Zero(params.back_stresses());
            Eigen::VectorXd Am = Eigen::VectorXd::Zero(params.back_stresses());
            while(abs(dDL) > 1e-15) {
                Vector6 Cij = st;
                D = params.sy0() + 3*G*DL;

                Vector6 dCijdDL = Vector6::Zero();
                double dDdDL = 3*G + dsydDL;

                if (params.isotropic_hardening()) {
                    R2 = (state.R() + params.b()*params.Q()*DL)/(1 + params.b()*DL);
                    sy = params.sy0() + R2;
                    D += R2;
                    dsydDL = params.b()*(params.Q() - state.R());
                }

                if (params.kinematic_hardening()) {
                    for (unsigned i=0; i != params.back_stresses(); ++i) {
                        double Cm = params.Cm(i);
                        double gm = params.gamma(i);

                        theta[i] = 1./(1+gm*DL);
                        Am[i] = gm/((1+gm*DL)*(1+gm*DL));

                        Cij -= theta[i]*state.back_stress_vector(i);
                        std::cout << Cij << std::endl;
                        D += DL*theta[i]*Cm;
                        std::cout << "Cij=" << Cij << " D=" << D << std::endl;
                        dCijdDL += Am[i]*state.back_stress_vector(i);
                        dDdDL += theta[i]*theta[i]*Cm;
                    }
                }

                nij = 1.5*Cij/D;
                if (phase_transformations) {
                    RA = params.R1() + params.R2()*sy/params.sy0_au();
                    D += 3*G*RA*DfM;
                    dDdDL += params.R2()/params.sy0_au()*DfM*dsydDL;
                }

                double f = 2./3*double_contract(nij, nij) - 1;
                Vector6 dndDL = -nij*dDdDL/D;
                double dfdDL = 4./3*double_contract(nij, dndDL);
                if (! phase_transformations) {
                    dDL = f/dfdDL;
                    DL -= dDL;
                }
            }

            state.ep_eff() += DL;
            state.R() = R2;
            stress_vec -= 2*G*DL*nij;
            D_alg = Del;
        }
    }
}