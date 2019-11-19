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

double ms_stress(const Eigen::Matrix<double, 6, 1>& stress, const TransformationMaterialParameters& params) {
    Eigen::Matrix<double, 6, 1> s_dev = deviator(stress);
    double m_stress = params.a1()*(stress[0] + stress[1] + stress[2]);   // Contribution from hydrostatic stress
    m_stress += params.a2()*von_Mises(stress);
    m_stress += params.a3()*vector_det(s_dev);
    return m_stress;
}

double ms_strain(double epl, const TransformationMaterialParameters& params) {
    return params.beta()*pow(1 - exp(-params.alpha()*epl), 4);
}

double transformation_function(const Eigen::Matrix<double, 6, 1>& stress, double epl, double T,
                               const TransformationMaterialParameters& params) {
    double a = exp(-params.k()*(params.Ms() + ms_stress(stress, params) +
                                ms_strain(epl, params) + params.Mss() - T));
    return 1 - a;
}

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
        double *rpl, double *ddsddt, double *drplde, double *drpldt, double *stran, double *dstran, double* time,
        double& dtime, double& temp, double& dtemp, double *predef, double *dpred, char *cmname, const int& ndi,
        const int& nshr, const int& ntens, const int& nstatv, const double* props, const int& nprops, double *coords,
        double* drot, double *pnewdt, double& celent, double* dfgrd0, double* dfgrd1, const int& noel, const int& npt,
        const int& layer, const int& kspt, const int& kstep, const int& kinc, short cmname_len) {

    using Matrix6x6 = Eigen::Matrix<double, 6, 6>;
    using Vector6 = Eigen::Matrix<double, 6, 1>;

    const TransformationMaterialParameters params(props);
    State state(statev, params.back_stresses());

    Eigen::Map<Vector6> stress_vec(stress);

    const Eigen::Map<Vector6> de(dstran);
    Eigen::Map<Matrix6x6> D_alg(ddsdde);

    // Elastic parameters
    double G = params.E()/2/(1 + params.v());
    double K = params.E()/3/(1 - 2*params.v());
    Matrix6x6 Del = 2*G*J + K*E3;

    // Yield stress at start of increment
    double sy = params.sy0M()*state.fM() + params.sy0A()*(1-state.fM()) + state.R();

    Vector6 sigma_t = stress_vec + Del*de;  // Trial stress
    Vector6 sij_t = deviator(sigma_t);

    Vector6 stilde = sij_t;
    if (params.kinematic_hardening()) {
        stilde -= state.total_back_stress();
    }
    bool plastic = params.plastic() && yield_function(stilde, sy) > 0;
    std::cout << temp << "  " << transformation_function(sigma_t, state.ep_eff(), temp, params) - state.fM() <<  std::endl;
    bool phase_transformations = transformation_function(sigma_t, state.ep_eff(), temp, params) - state.fM() >= 0;
    std::cout << "phase_transformations " << phase_transformations << std::endl;
    bool elastic = !plastic && !phase_transformations;
    if (elastic) {     // Use the trial stress as the stress and the elastic stiffness matrix as the tangent
        D_alg = Del;
        stress_vec = sigma_t;
    }
    else {  // Inelastic deformations
        // Increment in plastic strain and martensitic phase fraction
        Vector6 sigma_2 = sigma_t;

        double DL = 0;
        double DfM = 0;

        double f = 0;
        double h = 0;

        double dfdDL = 0;
        double dfdDfM = 0;
        double dhdDL = 0;
        double dhdDfM = 0;

        // Effective stress and its derivatives
        double s_eq_2 = sqrt(1.5*double_contract(stilde, stilde));
        double ds_eq_2_dDL = 0;
        double ds_eq_2_dfM = 0;

        double dMepdDL = 0;

        double dR2dDL = 0;
        double RA = 0;
        double s_eq_prime = s_eq_2;
        Matrix6x6 nnt = Matrix6x6::Zero();
        Matrix6x6 Aijkl = Matrix6x6::Zero();
        Vector6 nij2 = 1.5*stilde;
        if (s_eq_2 > 0) {
            nij2 /= s_eq_2;
        }

        double F = 0;
        Vector6 bij = Vector6::Zero();

        double R2 = 0;

        Vector6 dsij_prime_dDL = Vector6::Zero();
        double B = 1;
        double residual = 1e99;
        unsigned iter = 0;
        while (residual > 1e-15) {
            ++iter;
            sigma_2 = sigma_t;

            double dDL = 0;
            double dDfM = 0;

            B = 1 + 3*G*params.R2()*DfM/params.sy0A();
            s_eq_2 = (s_eq_prime - 3*G*(DL + params.R1()*DfM))/B;
            if (plastic) {
                dsij_prime_dDL = Vector6::Zero();
                ds_eq_2_dDL = -3*G;
                double sy0 = params.sy0M()*(state.fM() + DfM) + params.sy0A()*(1 - (state.fM() + DfM));
                R2 = (state.R() + params.b()*params.Q()*DL)/(1 + params.b()*DL);
                double sy_2 = sy0 + R2;
                Vector6 sij_prime = sij_t;
                double back_stress_correction = 0;
                for (unsigned i = 0; i != params.back_stresses(); ++i) {
                    double theta = 1./(1 + params.gamma(i)*DL);
                    sij_prime -= theta*state.back_stress_vector(i);
                    back_stress_correction += theta*params.Cm(i)*DL;
                    dsij_prime_dDL += params.gamma(i)*theta*theta*state.back_stress_vector(i);
                    ds_eq_2_dDL -= theta*theta*params.Cm(i);
                }
                s_eq_prime = sqrt(1.5*double_contract(sij_prime, sij_prime));
                s_eq_2 = (s_eq_prime - 3*G*(DL + params.R1()*DfM) - back_stress_correction)/B;
                nij2 = 1.5*sij_prime/s_eq_prime;
                ds_eq_2_dDL += double_contract(nij2, dsij_prime_dDL);
                ds_eq_2_dDL /= B;
                dR2dDL = params.b()/(1 + params.b()*DL)*(params.Q() - R2);
                dfdDL = ds_eq_2_dDL - dR2dDL;
                f = s_eq_2 - sy_2;
                sigma_2 -= 2*G*DL*nij2;
            }
            if (phase_transformations) {
                RA = params.R1() + params.R2()*s_eq_2/params.sy0A();
                sigma_2 -= (2*G*RA*nij2 + K*params.dV()/3*delta_ij)*DfM;
                h = transformation_function(sigma_2, state.ep_eff() + DL, temp, params) - (state.fM() + DfM);
                F = params.k()*exp(-params.k()*(params.Ms() + ms_stress(sigma_2, params)
                                                + ms_strain(state.ep_eff() + DL, params) + params.Mss() - temp));
                Vector6 s = deviator(sigma_2);
                std::cout << "s: " << s.transpose().format(CleanFmt) << std::endl;

                double J2 = 0.5*double_contract(s, s);
                std::cout << "F " << F << std::endl;
                bij = params.a1()*delta_ij;
                if (J2 > 0) {
                    bij += 1.5*params.a2()*s/sqrt(3*J2) + params.a3()*(contract(s, s) - 2./3*J2*delta_ij);
                }
                bij *= F;
                ds_eq_2_dfM = -3*G*RA/B;
                Vector6 dsijdDfM = -2*G*(RA + DfM*params.R2()/params.sy0A()*ds_eq_2_dfM)*nij2 - K/3*params.dV()*delta_ij;
                dhdDfM = double_contract(bij, dsijdDfM) - 1;
                std::cout << "dsijdDfM: " << dsijdDfM.transpose().format(CleanFmt) << std::endl;
                std::cout << "bij: " << bij.transpose().format(CleanFmt) << std::endl;
                std::cout << "dhdDfM: " << dhdDfM << std::endl;
            }

            nnt = nij2*nij2.transpose();
            Aijkl = J - 2./3*nnt;
            if (plastic && phase_transformations) {
                dMepdDL = params.beta()*params.n()*pow((1 - exp(-params.alpha()*(state.ep_eff() + DL))),
                                                       params.n() - 1)*params.alpha()*
                          exp(-params.alpha()*(state.ep_eff() + DL));
                dfdDfM = -3*G*RA/B - (params.sy0M() - params.sy0A());
                Vector6 dsigmaijdDL = -2*G*(1 + DfM*params.R2()/params.sy0A()*ds_eq_2_dDL)*nij2;
                dhdDL = double_contract(bij, dsigmaijdDL) + F*dMepdDL;
                double detJ = dfdDL*dhdDfM - dfdDfM*dhdDL;
                dDL = (dhdDfM*-f - dfdDfM*-h)/detJ;
                dDfM = (-dhdDL*-f + dfdDL*-h)/detJ;
                if (dDL + DL < 0) {
                    dDL = 0;
                    DL = 0;
                    dDfM = -h/dhdDfM;
                }
                if (dDfM + DfM < 0) {
                    dDfM = 0;
                    DfM = 0;
                    dDL = -f/dfdDL;
                }
            } else if (plastic) {
                dDL = -f/dfdDL;
            } else {  // Only phase transformations
                dDfM = -h/dhdDfM;
            }
            DL += dDL;
            DfM += dDfM;
            residual = abs(dDL) + abs(dDfM);
        }
        // Updating state variables
        state.ep_eff() += DL;
        state.fM() += DfM;
        state.R() = R2;
        stress_vec = sigma_2;

        if (params.kinematic_hardening()) {
            state.total_back_stress() = Vector6::Zero();
            for (unsigned i = 0; i != params.back_stresses(); ++i) {
                state.back_stress_vector(i) += 2./3*params.Cm(i)*DL*nij2;
                state.back_stress_vector(i) /= (1 + params.gamma(i)*DL);
                state.total_back_stress() += state.back_stress_vector(i);
            }
        }
        D_alg = Del;
        if (s_eq_prime > 0) {
            D_alg -= 6*G*G*(DL + RA*DfM)/s_eq_prime*Aijkl;
        }
        std::cout << "Aijkl" << std::endl << Aijkl.format(CleanFmt) << std::endl;
        double A = dR2dDL - F*dMepdDL*dfdDfM -  ds_eq_2_dDL;
        Vector6 Lekl = 2*G/A/B*nij2;

        if (DL > 0) {
            D_alg -= 2*G*(1 + DfM*params.R2()/params.sy0A()*ds_eq_2_dDL)*nij2*Lekl.transpose();
        }

        if (DfM > 0) {
            D_alg -= 4*G*G/B*DfM*params.R2()/params.sy0A()*nnt;
            Matrix6x6 Bijkl = I;
            Vector6 Fskl = bij;
            if (DL > 0) {
                Fskl += F*dMepdDL/A*dfdDfM*bij;
                Vector6 Lskl = 1/A*dfdDfM*bij;
                Vector6 Fekl = F*dMepdDL/A/B*2*G*nij2;
                D_alg -= 2*G*(RA + DfM*params.R2()/params.sy0A()*ds_eq_2_dfM)*nij2*Fekl.transpose()
                        - K/3*params.dV()*delta_ij*Fekl.transpose();
                Bijkl += 2*G*(1+DfM*params.R2()/params.sy0A()*ds_eq_2_dDL)*nij2*Lskl.transpose();
            }
            Bijkl += 2*G*(RA + DfM*params.R2()/params.sy0A()*ds_eq_2_dfM)*nij2*Fskl.transpose()
                    + K/3*params.dV()*delta_ij*Fskl.transpose();
            for (unsigned i = 3; i != 6; ++i) {
                for (unsigned j = 3; j != 6; ++j)
                    Bijkl(i, j) *= 2;
            }
            D_alg = Bijkl.inverse()*D_alg;
            std::cout << D_alg.format(CleanFmt) << std::endl;
        }
    }
}
