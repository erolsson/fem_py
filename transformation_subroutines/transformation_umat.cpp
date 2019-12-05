//
// Created by erolsson on 11/09/2019.
//

#include "transformation_umat.h"

#include <cmath>
#include <iostream>

#include "Eigen/Dense"

#include "simulation_parameters.h"
#include "stress_functions.h"

const double pi = 3.14159265359;

class State {
public:
    using Vector6 = Eigen::Matrix<double, 6, 1>;
    explicit State(double* data, unsigned back_stresses) : data_(data), back_stresses_(back_stresses) {}
    double& ep_eff() { return data_ [0]; }
    double& fM() { return  data_[1]; }
    double& fM_stress() { return data_[2]; }
    double& fM_strain() { return data_[3]; }
    double& fsb() { return data_[4]; }
    double& fsb0() { return data_[5]; }
    double& R() { return data_ [6]; }


    Eigen::Map<Vector6> back_stress_vector(unsigned n) {
        return Eigen::Map<Vector6>(&data_[7 + n*6]);
    }

    Eigen::Map<Vector6> total_back_stress() {
        return Eigen::Map<Vector6>(&data_[7 + back_stresses_*6]);
    }

private:
    double* data_;
    unsigned back_stresses_;
};

double yield_function(const Eigen::Matrix<double, 6, 1>& sigma, const Eigen::Matrix<double, 6, 1>& alpha,
        double sigma_y, const TransformationMaterialParameters& params) {
    Eigen::Matrix<double, 6, 1> stilde = deviator(sigma) - alpha;
    double I1 = sigma[0] + sigma[1] + sigma[2];
    return sqrt(3*double_contract(stilde, stilde)/2) + params.a()*I1 - sigma_y;
}

double stress_temperature_transformation(const Eigen::Matrix<double, 6, 1>& stress,
                                         const TransformationMaterialParameters& params, double T) {
    Eigen::Matrix<double, 6, 1> s_dev = deviator(stress);
    double m_stress = params.a1()*(stress[0] + stress[1] + stress[2]);   // Contribution from hydrostatic stress
    m_stress += params.a2()*von_Mises(stress);
    m_stress += params.a3()*vector_det(s_dev);
    return params.k()*(params.Ms() + m_stress + params.Mss() - T);
}

double stress_transformation_function(const Eigen::Matrix<double, 6, 1>& stress, double T,
                                      const TransformationMaterialParameters& params, double fM) {
    return 1 - exp(-stress_temperature_transformation(stress, params, T)) - fM;
}

double normal_pdf(double x) {
    return exp(-0.5*x*x)/sqrt(2*pi);
}

extern "C" void umat_(double *stress, double *statev, double *ddsdde, double *sse, double *spd, double *scd,
        double *rpl, double *ddsddt, double *drplde, double *drpldt, double *stran, double *dstran, double* time,
        double& dtime, double& temp, double& dtemp, double *predef, double *dpred, char *cmname, const int& ndi,
        const int& nshr, const int& ntens, const int& nstatv, const double* props, const int& nprops, double *coords,
        double* drot, double& pnewdt, double& celent, double* dfgrd0, double* dfgrd1, const int& noel, const int& npt,
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

    Vector6 stilde2 = sij_t;
    if (params.kinematic_hardening()) {
        stilde2 -= state.total_back_stress();
    }
    bool plastic = params.plastic() && yield_function(sigma_t, state.total_back_stress(), sy, params) > 0;
    std::cout << yield_function(sigma_t, state.total_back_stress(), sy, params) << std::endl;
    bool stress_transformations = stress_transformation_function(sigma_t, temp, params, state.fM()) >= 0;
    bool strain_transformations = params.beta() > 0 && plastic;
    bool elastic = !plastic && !stress_transformations;
    if (elastic) {     // Use the trial stress as the stress and the elastic stiffness matrix as the tangent
        D_alg = Del;
        stress_vec = sigma_t;
    }
    else {  // Inelastic deformations
        // Increment in plastic strain and martensitic phase fraction
        Vector6 sigma_2 = sigma_t;

        double DL = 0;
        double DfM_stress = 0;
        double DfM_strain = 0;

        double DfM = 0;

        double f = 0;
        double h_stress = 0;
        double h_strain = 0;

        double dfdDL = 0;
        double dfdDfM = 0;
        double dh_stressDL = 0;
        double dh_stressDfM = 0;

        double dh_straindDL = 0;
        double dh_strainDfM = 0;

        double fsb2 = state.fsb();

        double As = 0;
        double Bs = 0;
        double DSigma = 0;
        double Sigma = 0;
        double DvM_inv = 0;
        double P = 0;
        double pdf = 0;
        double dfsb2dDL = 0;
        double norm_drivning_force = 0;
        Vector6 dsvMdsij = Vector6::Zero();
        double s_vM_2 = 0;

        // Effective stress and its derivatives
        Vector6 s1 = deviator(static_cast<Vector6>(stress_vec));
        double s_vM_1 = sqrt(1.5*double_contract(s1, s1));
        double s_eq_2 = sqrt(1.5*double_contract(stilde2, stilde2));

        double I1_1 = stress_vec[0] + stress_vec[1] + stress_vec[2];

        double ds_eq_2_dDL = 0;
        double ds_eq_2_dfM = 0;

        double dR2dDL = 0;
        double RA = 0;
        double s_eq_prime = s_eq_2;
        Matrix6x6 nnt = Matrix6x6::Zero();
        Matrix6x6 Aijkl = Matrix6x6::Zero();
        Vector6 nij2 = 1.5*stilde2;
        if (s_eq_2 > 0) {
            nij2 /= s_eq_2;
        }

        Vector6 bij = Vector6::Zero();

        double R2 = 0;

        Vector6 dsij_prime_dDL = Vector6::Zero();
        double B = 1;
        double residual = 1e99;
        unsigned iter = 0;
        while (residual > 1e-15) {
            ++iter;
            double fM2 = state.fM() + DfM;
            sigma_2 = sigma_t;

            double dDL = 0;
            double dDfM_stress = 0;
            double dDfM_strain = 0;

            B = 1 + 3*G*params.R2()*DfM/params.sy0A();
            // s_eq_2 = (s_eq_prime - 3*G*(DL + params.R1()*DfM))/B;
            double I1_2 = sigma_t[0] + sigma_t[1] + sigma_t[2] - 3*K*params.dV()*DfM;

            // Calculates f and the derivative df/dDL,
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
                f = s_eq_2 + params.a()*I1_2 - sy_2;
            }
            RA = params.R1() + params.R2()*s_eq_2/params.sy0A();

            dfdDfM = -3*G*RA/B - params.a()*K*params.dV() - (params.sy0M() - params.sy0A());
            Vector6 dsijdDfM = -K*params.dV()*delta_ij;
            ds_eq_2_dfM = -3*G*RA/B;
            if ( s_eq_prime > 1e-12) {
                dsijdDfM -= 2*G*(RA + DfM*params.R2()/params.sy0A()*ds_eq_2_dfM)*nij2;
                sigma_2 -= 2*G*(DL + RA*DfM)*nij2 + K*params.dV()*delta_ij*DfM;
            }
            Vector6 dsijdDL = -2*G*(1 + DfM*params.R2()/params.sy0A()*ds_eq_2_dDL)*nij2;
            // Calculating the von Mises stress at step 2
            Vector6 s = deviator(sigma_2);
            double J2 = 0.5*double_contract(s, s);
            s_vM_2 = sqrt(3*J2);
            dsvMdsij = Vector6::Zero();
            if (J2 > 1e-12)
                dsvMdsij = 1.5*s/sqrt(3*J2);

            // h_strain and derivatives of h_strain
            if (plastic) {
                if (abs(s_vM_2 - s_vM_1) > 0) {
                    DvM_inv = 1/(s_vM_2 - s_vM_1);
                }
                Sigma = I1_2/s_vM_2;
                DSigma = (I1_2 - I1_1)*DvM_inv;

                double dSigmadDL = -Sigma/s_vM_2*double_contract(dsvMdsij, dsijdDL);
                double dSigmadDfM = 1/s_vM_2*(3*K*params.dV() + Sigma*double_contract(dsvMdsij, dsijdDfM));
                double dDSigmadDL = 0;
                double dDSigmadDfM = 0;
                    dDSigmadDL = -DSigma*DvM_inv*double_contract(dsvMdsij, dsijdDL);
                    dDSigmadDfM = DvM_inv*(3*K*params.dV() + DSigma*double_contract(dsvMdsij, dsijdDfM));
                fsb2 = 1 - (1 - state.fsb0())*exp(-params.alpha()*(state.ep_eff() + DL));
                dfsb2dDL = params.alpha()*(1 - state.fsb0())*exp(-params.alpha()*(state.ep_eff() + DL));

                double c = params.alpha()*params.beta()*(1-fsb2)*pow(fsb2, params.n()-1);
                double dcdDL = params.alpha()*params.beta()*((params.n()-1)*pow(fsb2, params.n()-2)
                        - params.n()*pow(fsb2, params.n()-1));

                double Gamma = params.g0() - params.g1()*temp/params.Ms() + params.g2()*Sigma;
                norm_drivning_force = (Gamma - params.g_mean())/params.g_std();
                P = 0.5*(1 + erf(norm_drivning_force));
                As = c*P;
                pdf = normal_pdf(norm_drivning_force)/params.g_std();
                Bs = params.g2()*params.beta()*pow(fsb2, params.n())*pdf*(DSigma > 0);
                double dAsddL = dcdDL*P + c*params.g2()*pdf*dSigmadDL;
                double dAsdfM = c*params.g2()*pdf*dSigmadDfM;

                double dBsdDL = params.g2()*params.beta()*pdf*pow(fsb2, params.n() - 1)*(params.n()*dfsb2dDL
                        - params.g2()*norm_drivning_force/params.g_std()*dSigmadDL)*(DSigma > 0);
                double dBsdfM = params.g2()*params.beta()*pow(fsb2,
                        params.n())*pdf*norm_drivning_force/params.g_std()*params.g2()*dSigmadDfM*(DSigma > 0);
                // std::cout << "As: " << As << " Bs: " << Bs << " P: " << P << " pdf: " << pdf << " Gamma: "
                //           << Gamma << " Sigma:" << Sigma << std::endl;
                h_strain = (1 - fM2)*(As*DL + Bs*DSigma) - DfM_strain;
                dh_straindDL = (1 - fM2)*(As + DL*dAsddL + Bs*dDSigmadDL + DSigma*dBsdDL);
                dh_strainDfM = -(As*DL + Bs*DSigma) + (1 - fM2)*(DL*dAsdfM + Bs*dDSigmadDfM + DSigma*dBsdfM);
                std::cout << "f: " << f << " dfdDL: " << dfdDL << " dfDfM: " << dfdDfM << std::endl;
                std::cout << "h_strain: " << h_strain << " dh_straindDL: " << dh_straindDL << " dh_strainDfM: "
                           << dh_strainDfM << std::endl;
            }

            if (stress_transformations) {
                h_stress = stress_transformation_function(sigma_2, temp, params, fM2);

                bij = params.a1()*delta_ij;

                if (J2 > 1e-12) {
                    bij += params.a2()*dsvMdsij+ params.a3()*(contract(s, s) - 2./3*J2*delta_ij);
                }
                double exp_fun = 1 - (h_stress + fM2);
                bij *= exp_fun*params.k();
                dh_stressDfM = double_contract(bij, dsijdDfM) - 1;
            }

            if ( !plastic ) {
                dDfM_stress = h_stress/dh_stressDfM;

            }

            else if ( !stress_transformations) {
                if ( !strain_transformations || h_strain == 0.) {
                    dDL = f/dfdDL;
                }
                else {
                    double det = dfdDL*(dh_strainDfM - 1) - dfdDfM*dh_straindDL;
                    dDL = ((dh_straindDL-1)*f - dfdDfM*h_strain)/det;
                    dDfM_strain = (-dh_straindDL*f + dfdDL*h_strain)/det;
                }
            }
            else {
                double a = dfdDL;
                double b = dfdDfM;
                double c = dh_stressDL;
                double d = dh_stressDfM;
                double e = dh_straindDL;
                double g = dh_strainDfM;
                double det = a*d- c*b;   // Pseudo determinant appearing in the inverse
                dDL = (d*f - b*h_stress)/det;
                dDfM_stress = -(c + d*e -c*g)/det*f + (a + b*e - a*g)*h_stress + h_strain;
                dDfM_strain = (d*e - c*g)/det*f - (b*e-a*g)/det*h_stress - h_strain;
            }

            DL -= dDL;
            DfM_stress -= dDfM_stress;
            DfM_strain -= dDfM_strain;
            DfM = DfM_stress + DfM_strain;
            residual = abs(dDL) + abs(dDfM_stress) + abs(dDfM_strain);

            std::cout << "DL: " << DL << " DfM_stress: " << DfM_stress << " DfM_strain: " << DfM_strain << std::endl;
            if (iter > 25) {
                pnewdt = 0.25;
                return;
            }
        }
        std::cout << "Converged in: " << iter << "iterations"  << std::endl;
        std::cout << "DL: " << DL << " DfM_stress: " << DfM_stress << " DfM_strain: " << DfM_strain
                  << " DfM: " << DfM << std::endl;
        // Updating state variables
        state.ep_eff() += DL;
        state.fM() += DfM;
        state.fM_stress() += DfM_stress;
        state.fM_strain() += DfM_strain;
        state.R() = R2;
        state.fsb() = fsb2;
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

        double A = dR2dDL -  ds_eq_2_dDL;

        Vector6 Lekl = Vector6::Zero();
        Vector6 Lskl = Vector6::Zero();
        Vector6 Fekl = Vector6::Zero();
        Vector6 Fskl = Vector6::Zero();


        if (DfM_stress > 0) {
            Fskl = bij;
            Lekl = 1./A/B*(2*G*nij2 + params.a()*K*delta_ij);
            Lskl = 1/A*dfdDfM*bij;
        }
        else {
            double B1 = (1 - state.fM())/(1 + As*DL + Bs*DSigma);
            double B2 = B1*(As + params.alpha()*params.beta()*pow(fsb2, params.n()-2)*(params.n()*(1-fsb2) - 1)*P
                  + pdf*params.n()*params.beta()*pow(fsb2, params.n() - 1)*(DSigma > 0)*dfsb2dDL);
            Lekl = (2*G*nij2 + params.a()*K*delta_ij)/(A*B - B*B2*dfdDfM);

            Vector6 Gkl = B1*Bs*(DvM_inv*delta_ij - DSigma*dsvMdsij +
                    (norm_drivning_force/params.g_std()*params.g2()*DSigma*(DSigma > 0)
                    *(1/s_vM_2*delta_ij - Sigma*dsvMdsij)));
            Lskl = Gkl*dfdDfM/(A + B2*dfdDfM);
            Fekl = B2*Lekl;
            Fskl = Gkl*(1 + dfdDfM/(A + B2*dfdDfM)*B2);
        }

        if (DL > 0) {
            D_alg -= 2*G*(1 + DfM*params.R2()/params.sy0A()*ds_eq_2_dDL)*nij2*Lekl.transpose();
        }

        if (DfM > 0) {
            D_alg -= 4*G*G/B*DfM*params.R2()/params.sy0A()*nnt;
            Matrix6x6 Bijkl = I;
            if (DL > 0) {
                D_alg -= 2*G*(RA + DfM*params.R2()/params.sy0A()*ds_eq_2_dfM)*nij2*Fekl.transpose()
                        - K*params.dV()*delta_ij*Fekl.transpose();
                Bijkl += 2*G*(1+DfM*params.R2()/params.sy0A()*ds_eq_2_dDL)*nij2*Lskl.transpose();
            }
            Bijkl += 2*G*(RA + DfM*params.R2()/params.sy0A()*ds_eq_2_dfM)*nij2*Fskl.transpose()
                    + K*params.dV()*delta_ij*Fskl.transpose();
            for (unsigned i = 3; i != 6; ++i) {
                for (unsigned j = 3; j != 6; ++j)
                    Bijkl(i, j) *= 2;
            }
            D_alg = Bijkl.inverse()*D_alg;
            std::cout << D_alg.format(CleanFmt) << std::endl;
        }
    }
}
