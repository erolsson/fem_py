//
// Created by erolsson on 11/09/2019.
//

#ifndef UTILITIES_H
#define UTILITIES_H

#include "Eigen/Dense"

using Eigen::Matrix;
// Unit forth order tensor in matrix form
const static Matrix<double, 6, 6> I((Matrix<double, 6, 6, Eigen::RowMajor>() <<
        1.,  0.,  0.,  0.,  0.,  0.,
        0.,  1.,  0.,  0.,  0.,  0.,
        0.,  0.,  1.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.5, 0.,  0.,
        0.,  0.,  0.,  0.,  0.5, 0.,
        0.,  0.,  0.,  0.,  0.,  0.5).finished());

//Deviatoric fourth order tensor in matrix form
const static Matrix<double, 6, 6> J((Matrix<double, 6, 6, Eigen::RowMajor>() <<
         2./3,  -1./3,  -1./3,  0.,  0.,  0.,
        -1./3.,  2./3,  -1./3,  0.,  0.,  0.,
        -1./3,  -1./3,   2./3., 0.,  0.,  0.,
         0.,     0.,     0.,    0.5, 0.,  0.,
         0.,     0.,     0.,    0.,  0.5, 0.,
         0.,     0.,     0.,    0.,  0.,  0.5).finished());

const static Matrix<double, 6, 6> E3((Matrix<double, 6, 6, Eigen::RowMajor>() <<
        1.,  1.,  1.,  0.,  0.,  0.,
        1.,  1.,  1.,  0.,  0.,  0.,
        1.,  1.,  1.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.).finished());

#endif //UTILITIES_H
