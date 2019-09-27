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

const static SimulationParameters* props;
/*
double yield_function(const Eigen::Matrix<double, 6, 1>& stilde, double sigma_y) {
    return double_contract(stilde, stilde) - sigma_y*sigma_y;
}

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
extern "C" void uexternaldb_(const int* lop, const int* lrestart, const double* time, const double* dtime,
                             const int* kstep, const int* kinc) {
    char out_dir_char[256];
    int out_dir_len;
    getoutdir_(out_dir_char, out_dir_len, 256);
    std::string out_dir(out_dir_char, out_dir_char+out_dir_len);

    char job_name_char[256];
    int job_name_len;
    getjobname_(job_name_char, job_name_len, 256);
    std::string job_name(job_name_char, job_name_char + job_name_len);
    std::string matierial_file_name = out_dir + "/" + job_name + ".par";
    std::fstream outfile(matierial_file_name);

    if (!outfile.good()) {
        matierial_file_name = out_dir + "/material_parameters.par";
        outfile = std::fstream(matierial_file_name);
        if (!outfile.good()) {
            std::cerr << "No material_parameters.par or " << job_name << ".par in the running directory" << std::endl;
            std::cerr << "Exiting!" << std::endl;
            std::abort();
        }
    }

    if (*lop == 0) {
        props = new SimulationParameters(matierial_file_name);
    }
    else if (*lop == 3) {
        delete props;
    }
}

extern "C" void umat_(double *stress, double *statev, double *ddsdde, double *sse, double *spd, double *scd,
        double *rpl, double *ddsddt, double *drplde, double *drpldt, double *stran, double *dstran, double *time,
        double *dtime, double *temp, double *dtemp, double *predef, double *dpred, char *cmname, int *ndi, int *nshr,
        int *ntens, int *nstatv, const double*, int *nprops, double *coords, double *drot, double *pnewdt,
        double *celent, double *dfgrd0, double *dfgrd1, int *noel, int *npt, int *layer, int *kspt, int *kstep,
        int *kinc, short cmname_len) {

    using Matrix6x6 = Eigen::Matrix<double, 6, 6, Eigen::RowMajor>;
    using Vector6 = Eigen::Matrix<double, 6, 1>;

    double G = props->E/2/(1+props->v);
    double K = props->E/3/(1-2*props->v);
    // Collecting state variables
    double R = statev[0];
    // double sy = sy0 + R;
    // Vector6 back_stress;
    // back_stress << statev[1], statev[2], statev[3], statev[4], statev[5], statev[6];
    Eigen::Map<Matrix6x6> D_alg(ddsdde);
    Matrix6x6 Del = 2*G*J + K*E3;
    D_alg = Del;

    Eigen::Map<Vector6> s1(stress);
    Vector6 de = Eigen::Map<Vector6>(dstran);
    s1 += Del*de;           // Trial stress

    // Vector6  stilde = deviator(st) - back_stress;
    // bool plastic = yield_function(stilde, sy) > 0;

}