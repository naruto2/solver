#ifndef CULU_H
#define CULU_H
#include <cusolverDn.h> // dense LAPACK

typedef struct {
  double *workspace;
  double *dB;
  int n;
  int nrhs;
  cusolverDnHandle_t handle;
  double *dA;
  int lda;
  int *pivot;
  int *devInfo;
} LUparameter;


vector<double> LUsolve(LUparameter LUp, const vector<double> &b);
void LUdestroy(LUparameter LUp);
LUparameter LUinit(matrix<double> AA);

#endif
