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
    return 0.;
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
    // std::cout << "de: " << de.transpose().format(CleanFmt) << std::endl;
    // std::cout << "sigma_0: " << stress_vec.transpose().format(CleanFmt) << std::endl;
    // std::cout << "sigma_t: " << sigma_t.transpose().format(CleanFmt) << std::endl;
    Vector6 sij_t = deviator(sigma_t);

    Vector6 stilde = sij_t;
    if (params.kinematic_hardening()) {
        stilde -= state.total_back_stress();
    }
    bool plastic = params.plastic() && yield_function(stilde, sy) > 0;
    bool phase_transformations = transformation_function(sigma_t, 0, temp, params) - state.fM() > 1e-12;
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
        //std::cout << "seq2: " << s_eq_2 << " sy:"  << sy << std::endl;

        double F = 0;
        Vector6 bij = Vector6::Zero();

        double R2 = 0;

        Vector6 dsij_prime_dDL = Vector6::Zero();
        double B = 1;
        double residual = 1e99;
        unsigned iter = 0;
        while(residual > 1e-15) {
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
            std::cout << "time:" << time[0] << std::endl;
            std::cout << "fM: " << state.fM() << std::endl;
            std::cout << "st: " << sigma_t.transpose().format(CleanFmt) << std::endl;
            if (phase_transformations) {
                std::cout << "s_eq_2: " << s_eq_2 << std::endl;
                RA = params.R1() + params.R2()*s_eq_2/params.sy0A();
                std::cout << "RA: " << std::endl;
                sigma_2 -= (2*G*RA*nij2 + K*params.dV()/3*delta_ij)*DfM;
                std::cout << "sigma_2: " << sigma_2.transpose().format(CleanFmt) << std::endl;
                Vector6 stemp = sigma_t - (2*G*RA*nij2 + K*params.dV()/3*delta_ij)*1e-4;
                std::cout << "sigma_temp: " << stemp.transpose().format(CleanFmt) << std::endl;
                std::cout << "num_d: " << ((stemp - sigma_2)/1e-4).transpose().format(CleanFmt) << std::endl;
                h = transformation_function(sigma_2, state.ep_eff() + DL, temp, params) - (state.fM() + DfM);
                F = params.k()*exp(-params.k()*(params.Ms() + ms_stress(sigma_2, params)
                                + ms_strain(state.ep_eff() + DL, params) + params.Mss() - temp));
                std::cout << "F: " << F << std::endl;
                Vector6 s = deviator(sigma_2);
                double J2 = 0.5*double_contract(s, s);
                bij = F*(params.a1()*delta_ij + 1.5*params.a2()*s/sqrt(3*J2)
                         + params.a3()*(contract(s, s) - 2./3*J2*delta_ij));
                std::cout << "b: " << bij.transpose().format(CleanFmt) << std::endl;
                ds_eq_2_dfM = -3*G*RA/B;
                Vector6 dsijdDfM = -2*G*(RA + DfM*params.R2()/params.sy0A()*ds_eq_2_dfM)*nij2 - K/3*params.dV()*delta_ij;
                std::cout << "dsijdDfM: " << dsijdDfM.transpose().format(CleanFmt) << std::endl;
                dhdDfM = double_contract(bij, dsijdDfM) - 1;
                std::cout << "dhdDfM: " << dhdDfM << std::endl;
            }

            nnt = nij2*nij2.transpose();
            Aijkl = J - 2./3*nnt;
            if (plastic && phase_transformations) {
                dfdDfM = -3*G*params.R1()/B - (params.sy0M() - params.sy0A());
                Vector6 dsigmaijdDL = -3*G*double_contract(Aijkl, dsij_prime_dDL)*0
                        - 2*G*(1 + 0*DfM*params.R2()/params.sy0A()*ds_eq_2_dDL)*nij2;
                dhdDL = double_contract(bij, dsigmaijdDL) + F*dMepdDL;
                double detJ = dfdDL*dhdDfM - dfdDfM*dhdDL;
                dDL = (dhdDfM*-f - dfdDfM*-h)/detJ;
                dDfM = (-dhdDL*-f + dfdDL*-h)/detJ;
                std::cout << DfM*params.R2()/params.sy0A()*ds_eq_2_dDL<< std::endl;
                std::cout << "bij: " << bij.transpose().format(CleanFmt) << std::endl;
                std::cout << "ds_eq_2_dDL: " << ds_eq_2_dDL << std::endl;
                std::cout <<  " a2: "
                          << (2*G*(1 + DfM*params.R2()/params.sy0A()*ds_eq_2_dDL)*nij2).transpose().format(CleanFmt)
                          << std::endl;
                std::cout << "J = " << std::endl << dfdDL << ", " << dfdDfM << std::endl << dhdDL << ", " << dhdDfM << std::endl;
            }
            else if (plastic) {
                dDL = -f/dfdDL;
            }
            else {  // Only phase transformations
                dDfM = -h/dhdDfM;
            }
            std::cout << "DL:" << DL << "  DfM:" << DfM << std::endl;
            std::cout << "dDL:" << dDL << "  dDfM:" << dDfM << std::endl;
            std::cout << "f: " << f << "  h: " << h << std::endl;
            DL += dDL;
            DfM += dDfM;
            residual = abs(dDL) + abs(dDfM);
        }
        std::cout << "Converged in " << iter << " iterations" << std::endl;
        // Updating state variables
        state.ep_eff() += DL;
        state.fM() += DfM;
        state.R() = R2;
        stress_vec = sigma_2;

        if (params.kinematic_hardening()) {
            state.total_back_stress() = Vector6::Zero();
            for (unsigned i = 0; i != params.back_stresses(); ++i) {
                state.back_stress_vector(i) += 2./3*params.Cm(i)*DL*nij2;
                state.back_stress_vector(i) /= (1+params.gamma(i)*DL);
                state.total_back_stress() += state.back_stress_vector(i);
            }
        }
        D_alg = Del;
        D_alg -=  6*G*G*(DL + RA*DfM)/s_eq_prime*Aijkl;
        double A = dR2dDL - F*dMepdDL*dfdDfM -  ds_eq_2_dDL;

        if (DL > 0) {
            D_alg -= 4*G*G/A/B*nnt
                    + 6*G*G*(DL + RA*DfM)/s_eq_prime*(1./A/B*double_contract(Aijkl, dsij_prime_dDL)*nij2.transpose());
            D_alg -= 4*G*G*(DfM*params.R2()/B/params.sy0A()*(ds_eq_2_dDL + ds_eq_2_dfM*F*dMepdDL))*nnt;
        }

        if (DfM > 0) {
            D_alg -= 4*G*G/B*DfM*params.R2()/params.sy0A()*nnt;
            Matrix6x6 Bijkl = I + K/3*params.dV()*delta_ij*bij.transpose();
            double B1 = RA + DfM*params.R2()/params.sy0A()*ds_eq_2_dfM;

            if (DL > 0) {
                B1 += 1/A*dfdDfM*(1+F*dMepdDL + dfdDfM*params.R2()/params.sy0A()*(ds_eq_2_dDL + ds_eq_2_dfM*F*dMepdDL));
                Bijkl += 3*G*(DL + RA*DfM)/s_eq_prime/A*dfdDfM*double_contract(Aijkl, dsij_prime_dDL)*bij.transpose();
            }

            Bijkl += 2*G*B1*nij2*bij.transpose();
            for (unsigned i = 3; i != 6; ++i) {
                Bijkl(i, i) *= 2;
            }
            D_alg = Bijkl.inverse()*D_alg;
        }
    }
}