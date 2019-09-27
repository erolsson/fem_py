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
    double& R() { return data_ [1]; }

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

    if (elastic) {     // Use the trial stress as the stress and the elastic stiffness matrix as the tangent
        stress_vec = st;
        D_alg = Del;
    }
    else {
        if (plastic) {
            // Calculating the increment in effective plastic strain DL
            double DL = 0;
            double dDL = 1e99;
            double sy2 = params.sy0();
            while(abs(dDL) < 1e-15) {

                double R = (state.R() + params.b()*params.Q()*DL)/(1 + params.b()*DL);
                if (params.isotropic_hardening()) {
                    sy2 += R;
                }

            }
        }
    }
}