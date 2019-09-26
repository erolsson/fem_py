//
// Created by erolsson on 11/09/2019.
//

#ifndef FEM_PY_TRANSFORMATION_UMAT_H
#define FEM_PY_TRANSFORMATION_UMAT_H

extern "C" void getoutdir_(char* outdir, int&, int);

extern "C" void getjobname_(char* jobname, int&, int);

extern "C" void umat_(double *stress, double *statev, double *ddsdde, double *sse, double *spd, double *scd,
                      double *rpl, double *ddsddt, double *drplde, double *drpldt, double *stran, double *dstran,
                      double *time, double *dtime, double *temp, double *dtemp, double *predef, double *dpred,
                      char *cmname, int *ndi, int *nshr, int *ntens, int *nstatv, const double *props, int *nprops,
                      double *coords, double *drot, double *pnewdt, double *celent, double *dfgrd0, double *dfgrd1,
                      int *noel, int *npt, int *layer, int *kspt, int *kstep, int *kinc, short cmname_len);

#endif //FEM_PY_TRANSFORMATION_UMAT_H
