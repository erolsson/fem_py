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
    explicit State(double* data) : data_(data) {}
    double& ep_eff() { return data_ [0]; }
    double& fM() { return  data_[1]; }
    double& R() { return data_ [2]; }

private:
    double* data_;
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
    State state(statev);

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
    bool elastic = !plastic && !phase_transformations;
    stress_vec = st;
    if (elastic) {     // Use the trial stress as the stress and the elastic stiffness matrix as the tangent
        D_alg = Del;
    }
    else {
        double DfM = 0;
        Vector6 nij;
        double RA = 0;
        if (plastic) {
            std::cout << "Plastic iteration" << std::endl;
            // Calculating the increment in effective plastic strain DL
            double DL = 0;
            double dDL = 1e99;

            double dsydDL = 0;

            double D = 0;
            double sy2 = sy;
            double R2 = 0;
            while(abs(dDL) < 1e-15) {
                double dDdDL = 3*G;
                D += 3*G*DL;
                if (params.isotropic_hardening()) {
                    R2 = (state.R() + params.b()*params.Q()*DL)/(1 + params.b()*DL);
                    sy2 = params.sy0() + R2;
                    D = sy2;
                    dsydDL = params.b()*(params.Q() - state.R());
                    dDdDL += (1 + params.R2()/params.sy0_au()*DfM)*dsydDL;
                }
                if (phase_transformations) {
                    RA = params.R1() + params.R2()*sy2/params.sy0_au();
                    D += 3*G*RA*DfM;
                }

                nij = 1.5*st/D;
                double f = 2./3*double_contract(nij, nij) - 1;

                Vector6 dnDL = -nij*dDdDL/D;
                double dfdDL = 4./3*double_contract(nij, dnDL);
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