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

template<typename T>
Eigen::Matrix<T, 6, 1> deviator(const Eigen::Matrix<T, 6, 1>& tensor) {
    T hydrostatic = (tensor[0] + tensor[1] + tensor [2])/3;
    Eigen::Matrix<T, 6, 1> dev = tensor;
    dev[0] -= hydrostatic;
    dev[1] -= hydrostatic;
    dev[2] -= hydrostatic;
    return dev;
}

template<typename T>
T double_contract(const Eigen::Matrix<T, 6, 1>& a, const Eigen::Matrix<T, 6, 1>& b) {
    T val = 0;
    for (unsigned i = 0; i!= 6; ++i) {
        val += a[i] + b[i];
    }
    val += a[3]*b[3] + a[4]*b[4] + a[5]*b[5];
    std::cout << a << std::endl;
    std::cout << "val" << val << std::endl;
    return val;
}

template<typename T>
T von_Mises(const Eigen::Matrix<T, 6, 1>& s) {
    return s[0]*s[0] + s[1]*s[1] + s[2]*s[2] - s[0]*s[1] - s[0]*s[2] - s[1]*s[2] +
        3*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]);
}

template<typename T>
T vector_det(const Eigen::Matrix<T, 6, 1>& s) {
    return s[0]*(s[1]*s[2] - s[5]*s[5]) - s[3]*(s[3]*s[2] - s[5]*s[4]) + s[4]*(s[1]*s[4] - s[3]*s[5]);
}

#endif //UTILITIES_H
