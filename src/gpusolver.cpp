#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cassert>
#include <vector>
#include <map>
#include "matrix.hpp"
using namespace std;


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
  int worksize;
} LUparameter;


template<typename T>
inline size_t bytesof(unsigned int s) { return s * sizeof(T); }

template<typename T>
T* allocate(unsigned int size) {
  T* result = nullptr;
  cudaError_t status = cudaMalloc(&result, bytesof<T>(size));
  assert( status == cudaSuccess );
  return result;
}


static vector<double> LUsolve(LUparameter&LUp, const vector<double> &b)
{
  static vector<double> x;
  x.resize(LUp.n);

  assert( LUp.handle != NULL );

  cusolverStatus_t status;
  
  cudaMemcpy(LUp.dB, &b[0], bytesof<double>(LUp.n*LUp.nrhs),
	     cudaMemcpyHostToDevice);

  status = cusolverDnDgetrs( LUp.handle, CUBLAS_OP_T, LUp.n, LUp.nrhs,
			     LUp.dA, LUp.lda, LUp.pivot, LUp.dB, LUp.n,
			     LUp.devInfo);
  assert( status == CUSOLVER_STATUS_SUCCESS );

  double *X = (double*)malloc(sizeof(double)*LUp.n);
  cudaMemcpy(X, LUp.dB, bytesof<double>(LUp.n*LUp.nrhs),
	     cudaMemcpyDeviceToHost);
  for(int i=0; i<LUp.n; i++) x[i] = X[i];
  free(X);
  return x;
}


static LUparameter &LUinit(matrix<double>& AA)
{
  static LUparameter LUp;
  cusolverStatus_t status;
  int n = LUp.n = AA.size();
  
  double *A = (double*)calloc(n*n,sizeof(double));
  assert( A != NULL );

  for (int i=0; i<n; i++ ) {
    auto AAi = AA[i];
    auto j = AAi.begin();
    for ( ; j != AAi.end(); j++ )
      A[i*n+j->first] = AA[i][j->first];
  }

  LUp.nrhs    = 1;
  LUp.lda     = n; 
  LUp.dA      = allocate<double>(n*n);
  LUp.dB      = allocate<double>(n*LUp.nrhs);
  LUp.pivot   = allocate<int>(n);
  LUp.devInfo = allocate<int>(1);
  
  status = cusolverDnCreate(&LUp.handle);
  assert( status == CUSOLVER_STATUS_SUCCESS );

  cudaMemcpy(LUp.dA, A, bytesof<double>(n*n), cudaMemcpyHostToDevice);
  //free(A);  A = NULL;
  
  status = cusolverDnDgetrf_bufferSize( LUp.handle, n, n, LUp.dA,
					LUp.lda, &LUp.worksize);
  assert( status == CUSOLVER_STATUS_SUCCESS );

  LUp.workspace = allocate<double>(LUp.worksize);
  
  status = cusolverDnDgetrf( LUp.handle, n, n, LUp.dA, LUp.lda,
			     LUp.workspace, LUp.pivot, LUp.devInfo);
  assert( status == CUSOLVER_STATUS_SUCCESS );

  return LUp;
}


static void LUdestroy(LUparameter &LUp)
{
  cudaFree(LUp.dA);              LUp.dA        = NULL;
  cudaFree(LUp.dB);              LUp.dB        = NULL;
  cudaFree(LUp.pivot);           LUp.pivot     = NULL;
  cudaFree(LUp.devInfo);         LUp.devInfo   = NULL;
  cudaFree(LUp.workspace);       LUp.workspace = NULL;
  cusolverDnDestroy(LUp.handle);
  cudaDeviceReset();
}


int gpusolver(matrix<double> &A, vector<double> &x,const vector<double> &b)
{
  LUparameter LUp = LUinit(A);
  x = LUsolve(LUp,b);
  LUdestroy(LUp);
  return 0;
}
