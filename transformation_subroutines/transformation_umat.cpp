//
// Created by erolsson on 11/09/2019.
//

#include "transformation_umat.h"

#include <iostream>

#include "Eigen/Dense"

#include "utilities.h"

extern "C" void umat_(double *stress, double *statev, double *ddsdde, double *sse, double *spd, double *scd,
        double *rpl, double *ddsddt, double *drplde, double *drpldt, double *stran, double *dstran, double *time,
        double *dtime, double *temp, double *dtemp, double *predef, double *dpred, char *cmname, int *ndi, int *nshr,
        int *ntens, int *nstatv, const double *props, int *nprops, double *coords, double *drot, double *pnewdt,
        double *celent, double *dfgrd0, double *dfgrd1, int *noel, int *npt, int *layer, int *kspt, int *kstep,
        int *kinc, short cmname_len)
{
    typedef Eigen::Matrix<double, 6, 6, Eigen::RowMajor> Matrix6x6;
    typedef Eigen::Matrix<double, 6, 1> Vector6;
    double E = props[0];
    double v = props[1];

    double G = E/2/(1+v);
    double K = E/3/(1-2*v);

    Eigen::Map<Matrix6x6> D_alg(ddsdde);
    Matrix6x6 Del = 2*G*J + K*E3;
    std::cout << 2*G*J << std::endl;
    D_alg = Del;
    std::cout << D_alg << std::endl;

    Eigen::Map<Vector6> s2 (stress);
    Vector6 de = Eigen::Map<Vector6>(dstran);
    s2 += Del*de;
}