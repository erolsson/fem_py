//
// Created by erolsson on 11/09/2019.
//

#include "transformation_umat.h"

#include <iostream>

#include "Eigen/Dense"

#include "simulation_parameters.h"
#include "stress_functions.h"

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
        double *dtime, double *temp, double *dtemp, double *predef, double *dpred, char *cmname, const int *ndi,
        const int *nshr, const int *ntens, const int *nstatv, const double* props, const int *nprops, double *coords,
        double *drot, double *pnewdt, double *celent, double *dfgrd0, double *dfgrd1, const int *noel, const int *npt,
        const int *layer, const int *kspt, const int *kstep, const int *kinc, short cmname_len) {

    using Matrix6x6 = Eigen::Matrix<double, 6, 6, Eigen::RowMajor>;
    using Vector6 = Eigen::Matrix<double, 6, 1>;

    const TransformationMaterialParameters params(props);
    double G = params.E()/2/(1+params.v());
    double K = params.E()/3/(1-2*params.v());

    // Collecting state variables
    State state(statev, params.back_stresses());

    double sy = params.sy0M()*state.fM() + params.sy0A()*(1-state.fM()) + state.R();
    // Vector6 back_stress;

    Eigen::Map<Matrix6x6> D_alg(ddsdde);
    Matrix6x6 Del = 2*G*J + K*E3;

    Eigen::Map<Vector6> stress_vec(stress);
    const Eigen::Map<Vector6> de(dstran);

    Vector6  st = stress_vec + Del*de;  // Trial stress

    Vector6 sij_t = deviator(st);
    Vector6 stilde = sij_t;
    if (params.kinematic_hardening()) {
        stilde -= state.total_back_stress();
    }
    bool plastic = params.plastic() && yield_function(stilde, sy) > 0;
    bool phase_transformations = params.strain_transformation || params.stress_transformation;
    phase_transformations = false;
    bool elastic = !plastic && !phase_transformations;
    stress_vec = st;
    if (elastic) {     // Use the trial stress as the stress and the elastic stiffness matrix as the tangent
        D_alg = Del;
    }
    else {  // Inelastic deformations
        double DfM = 0;
        double DL = 0;
        double RA = 0;
        double s_eq_prime = 0;

        Vector6 nij2;
        Vector6 st_dev = deviator(st);

        double residual = 1e99;

        double ds_eq_2_dDL = 0;
        double dsydDL = 0;

        double R2 = 0;

        Vector6 dsijdDL = Vector6::Zero();
        double B = 1;

        while(residual > 1e-15) {
            dsijdDL = Vector6::Zero();
            ds_eq_2_dDL = -3*G;
            double sy0 = params.sy0M()*(state.fM() + DfM) + params.sy0A()*(1-(state.fM() + DfM));
            R2 = (state.R() + params.b()*params.Q()*DL)/(1+params.b()*DL);
            double sy_2 = sy0 + R2;
            Vector6 sij_prime = sij_t;
            double back_stress_correction = 0;
            for (unsigned i = 0; i!= params.back_stresses(); ++i) {
                std::cout << params.Cm(i) << ", " << params.gamma(i) << ", " << state.back_stress_vector(i)[2] << std::endl;
                double theta = 1./(1+params.gamma(i)*DL);
                sij_prime -= theta*state.back_stress_vector(i);
                back_stress_correction += theta*params.Cm(i)*DL;
                dsijdDL += params.gamma(i)*theta*theta*state.back_stress_vector(i);
                ds_eq_2_dDL -= theta*theta*params.Cm(i);
            }
            s_eq_prime = sqrt(1.5*double_contract(sij_prime, sij_prime));
            B += 3*G*params.R2()*DfM/params.sy0A();
            double s_eq_2 = (s_eq_prime - 3*G*(DL+params.R1()*DfM) - back_stress_correction)/B;
            nij2 = 1.5*sij_prime/s_eq_prime;
            ds_eq_2_dDL += double_contract(nij2, dsijdDL);
            double dfdDL = ds_eq_2_dDL - params.b()/(1+params.b()*DL)*(params.Q() - R2);
            double f = s_eq_2 - sy_2;
            if (! phase_transformations) {
                double dDL = f/dfdDL;
                DL -= dDL;
                residual = abs(dDL);
            }
            std::cout << "DL=" << DL << std::endl;
        }

        // Updating state variables
        state.ep_eff() += DL;
        state.R() = R2;
        stress_vec -= 2*G*DL*nij2;

        if (params.kinematic_hardening()) {
            state.total_back_stress() = Vector6::Zero();
            for (unsigned i = 0; i != params.back_stresses(); ++i) {
                state.back_stress_vector(i) += 2./3*params.Cm(i)*DL*nij2;
                state.back_stress_vector(i) /= (1+params.gamma(i)*DL);
                std::cout << "da=" << 2./3*params.Cm(i)*DL*nij2[2] << " a=" << state.back_stress_vector(i)[2] << std::endl;
                state.total_back_stress() += state.back_stress_vector(i);
            }
        }

        if (DL > 0) {
            double A = - ds_eq_2_dDL;
            Matrix6x6 A_ijkl = J - 2./3*nij2*nij2.transpose();
            D_alg = Del - 4*G*G*nij2*nij2.transpose()/A - 6*G*G*DL/s_eq_prime*(A_ijkl -
                    1./A*double_contract(A_ijkl, dsijdDL)*nij2.transpose());
        }
    }
}