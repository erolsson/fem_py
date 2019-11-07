//
// Created by erolsson on 25/09/2019.
//

#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

#include "Eigen/Dense"

class TransformationMaterialParameters {
public:
    explicit TransformationMaterialParameters(const double* data):
            data_(data), back_stresses_(static_cast<int>(data[6])) {}

    // Elastic parameters
    [[nodiscard]] const double& E() const { return data_[0]; }         // Young's modulus
    [[nodiscard]] const double& v() const { return data_[1]; }         // Poisson's ratio

    // Initial yield stress of Martensite and Austenite
    [[nodiscard]] const double& sy0M() const { return data_[2]; } ;     // Initial yield stress
    [[nodiscard]] const double& sy0A() const { return data_[3]; } ;     // Initial yield stress

    // Parameters for isostropic hardnening
    [[nodiscard]] const double& Q() const { return data_[4]; } ;           // Maximum size of yield surface
    [[nodiscard]] const double& b() const { return data_[5]; } ;           // Growth rate of yield surface

    // Parameters for kinematic hardnening
    [[nodiscard]] unsigned back_stresses() const { return back_stresses_; }
    [[nodiscard]] const double& Cm(unsigned n) const { return data_[7 + 2*n]; }
    [[nodiscard]] const double& gamma(unsigned n) const { return data_[7 + 2*n+1]; }

    [[nodiscard]] const double& R1() const { return data_[7 + 2*back_stresses_]; } ;
    [[nodiscard]] const double& R2() const { return data_[8 + 2*back_stresses_]; } ;

    // Martensite start temperature
    [[nodiscard]] const double& Ms() const { return data_[9 + 2*back_stresses_]; } ;

    // Stabilisation temperature due to tempering
    [[nodiscard]] const double& Mss() const { return data_[10 + 2*back_stresses_]; } ;

    // Koistinen-Marburger parameter
    [[nodiscard]] const double& k() const { return data_[11 + 2*back_stresses_]; } ;

    // Material parameters for stress induced phase transformations
    [[nodiscard]] const double& a1() const { return data_[12 + 2*back_stresses_]; } ;
    [[nodiscard]] const double& a2() const { return data_[13 + 2*back_stresses_]; } ;
    [[nodiscard]] const double& a3() const { return data_[14 + 2*back_stresses_]; } ;

    [[nodiscard]] bool plastic() const { return sy0M() > 0 && sy0A() > 0; }
    [[nodiscard]] bool kinematic_hardening() const { return back_stresses_ > 0; }
    bool stress_transformation = false;
    bool strain_transformation = false;
private:
    const double* data_;
    unsigned back_stresses_ = 0;
};


#endif // SIMULATION_PARAMETERS_H
