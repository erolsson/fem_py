//
// Created by erolsson on 11/09/2019.
//

#include <iostream>

#include "transformation_umat.h"

int main() {
    double stress[6] = {0, 0, 0, 0, 0, 0};

    std::cout << stress[0] << ", " << stress[1] << ", " << stress[2] << ", "
              << stress[3] << ", " << stress[4] << ", " << stress[5] << std::endl;

    double ddsdde[36] = {0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0,};

    double dstran[6] = {0, 0, 1., 0, 0, 0};
    double props[2] = {200e3, 0.3};

    umat_(stress,     // stress
          0,          // statev
          ddsdde,     // ddsdde
          0,          // sse
          0,          // spd
          0,          // scd
          0,          // rpl
          0,          // ddsddt
          0,          // drplde
          0,          // drpldt
          0,          // stran
          dstran,     // dstran
          0,          // time
          0,          // dtime
          0,          // temp
          0,          // dtemp
          0,          // predef
          0,          // dpred
          0,          // cmname
          0,          // ndi
          0,          // nshr
          0,          // ntens
          0,          // nstatv
          props,      // props
          0,          // nprops,
          0,          // coords,
          0,          // drot,
          0,          // pnewdt,
          0,          // celent,
          0,          // dfgrd0,
          0,          // dfgrd1,
          0,          // noel,
          0,          // npt,
          0,          // layer,
          0,          // kspt,
          0,          // int *kstep,
          0,          // int *kinc,
          0      //short cmname_len
          );

    std::cout << stress[0] << ", " << stress[1] << ", " << stress[2] << ", "
              << stress[3] << ", " << stress[4] << ", " << stress[5] << std::endl;
    return 0;
}
