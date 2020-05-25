#ifndef SOLVER_H
#define SOLVER_H

#ifndef noGPU
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>
#include <thrust/device_vector.h>
using namespace thrust;
#endif

#include <iostream>
#include <cmath>
#include "solver/sparse.h"
#define SOLVERROR_NONE      0
#define SOLVERROR_MAXIT     1
#define SOLVERROR_BREAKDOWN 2
using namespace std;
using namespace sparse;
int progress(string str, long i, double res);
#include "solver/crs.h"
#include "solver/operator.h"


#include "solver/bicgstab.h"
#include "solver/ConjugateGradient.h"


#include "solver/mmio.h"
#include "solver/matrixmarket.h"

template < class Matrix >
int
symmetric(Matrix& A)
{
  return 0;
}

template < class Matrix, class Vector >
double
resi(Matrix& A, Vector& x, const Vector& b)
{
  Vector r = y_ax(-1.0*(A*x), 1.0, b);
  return nrm2(r)/nrm2(b);
}

template < class Matrix, class Vector >
int
solver(Matrix& A, Vector& x, const Vector&b)
{
  solveGMRES=0;
  M_init(A);
  Vector y = x; int ret = 1;
  for(int k=0; k< 16; k++){
    ret = QMR(A,x,b);
    if ( ret == 0 ) goto end;
    CRSinit(A);
    if (resi(A,x,b)>resi(A,y,b)) x=y; else y=x;
    CRSdestory(A);

    if(symmetric(A)) {
      ret = ConjugateGradient(A,x,b);
      if ( ret == 0 ) goto end;
      CRSinit(A);
      if (resi(A,x,b)>resi(A,y,b)) x=y; else y=x;
      CRSdestory(A);
    }

    ret = BiCG(A,x,b); 
    if ( ret == 0 ) goto end;
    CRSinit(A);
    if (resi(A,x,b)>resi(A,y,b)) x=y; else y=x;
    CRSdestory(A);

    ret = BiCGSTAB(A,x,b);
    if ( ret == 0 ) goto end; 
    CRSinit(A);
    if (resi(A,x,b)>resi(A,y,b)) x=y; else y=x;
    CRSdestory(A);
    
    ret = CGS(A,x,b); 
    if ( ret == 0 ) goto end;
    CRSinit(A);
    if (resi(A,x,b)>resi(A,y,b)) x=y; else y=x;
    CRSdestory(A);

    if ( k>0 ) solveGMRES=1;
    ret = GMRES(A,x,b); 
    if ( ret == 0 ) goto end;
    CRSinit(A);
    if (resi(A,x,b)>resi(A,y,b)) x=y; else y=x;
    CRSdestory(A);
  }
 end:
  M_destory(A);
  return ret;
}

int progress(string str, long i, double res)
{
  //cout<<str<<" i= "<<i<<" res= "<<res<<endl;
  if ( str == "QMR" ) {
    if ( res != res ) return 3;
  }
  if ( str == "BiCG" ) {
    if ( res > 10000.0 ) return 2; 
    if ( res != res ) return 3;
  }
  if ( str == "BiCGSTAB" ) {
    if ( res > 10000.0 ) return 2; 
    if ( res != res ) return 3;
  }
  if ( str == "CGS" ) {
    if ( res > 10000.0 ) return 2;
    if ( res != res ) return 3;
  }
  if ( str == "GMRES" ) {
    if ( res > 10000.0 ) return 2;
    if ( res != res ) return 3;
  }
  return 0;
}

#endif
