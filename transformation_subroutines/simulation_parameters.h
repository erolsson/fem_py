//
// Created by erolsson on 25/09/2019.
//

#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

#include "Eigen/Dense"

class TransformationMaterialParameters {
private:
    const double* data_;
    unsigned back_stresses_ = 0;
public:
    explicit TransformationMaterialParameters(const double* data):
            data_(data), back_stresses_(static_cast<int>(data[5])) {}

    // Elastic parameters
    [[nodiscard]] const double& E() const { return data_[0]; }         // Young's modulus
    [[nodiscard]] const double& v() const { return data_[1]; }         // Poisson's ratio

    //Plasticity
    [[nodiscard]] const double& sy0() const { return data_[2]; } ;     // Initial yield stress

    // Parameters for isostropic hardnening
    [[nodiscard]] const double& Q() const { return data_[3]; } ;           // Maximum size of yield surface
    [[nodiscard]] const double& b() const { return data_[4]; } ;           // Growth rate of yield surface
    [[nodiscard]] unsigned back_stresses() const { return back_stresses_; }
    [[nodiscard]] const double& Cm(unsigned n) const { return data_[6 + 2*n]; }
    [[nodiscard]] const double& gamma(unsigned n) const { return data_[6 + 2*n+1]; }

    [[nodiscard]] const double& R1() const { return data_[6 + 2*back_stresses_]; } ;
    [[nodiscard]] const double& R2() const { return data_[7 + 2*back_stresses_]; } ;
    [[nodiscard]] const double& sy0_au() const { return data_[8 + 2*back_stresses_]; } ;

    [[nodiscard]] bool plastic() const { return sy0() > 0; }
    [[nodiscard]] bool isotropic_hardening() const { return Q() > 0; }
    [[nodiscard]] bool kinematic_hardening() const { return back_stresses_ > 0; }
    bool stress_transformation = false;
    bool strain_transformation = false;

};


#endif // SIMULATION_PARAMETERS_H
