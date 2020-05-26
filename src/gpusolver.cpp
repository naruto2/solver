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


static vector<double> LUsolve(LUparameter LUp, const vector<double> &b)
{
  static vector<double> x;
  x.resize(LUp.n);

  if ( LUp.handle == NULL ) {
    for(int i=0; i<LUp.n; i++) x[i] = b[i];
    return x;
  }

  cusolverStatus_t status;
  
  cudaMemcpy(LUp.dB, &b[0], bytesof<double>(LUp.n*LUp.nrhs),
	     cudaMemcpyHostToDevice);
  int ldb = LUp.n;

  // AX = B を解く (解XはBをoverrideする)
  status = cusolverDnDgetrs(
             LUp.handle,
             CUBLAS_OP_T,
             LUp.n,     // 行(=列)
             LUp.nrhs,  // 問題数
             LUp.dA,    // A
             LUp.lda,   // Aのヨコハバ
             LUp.pivot, // LU分解で得られたピボット
             LUp.dB,    // B
             ldb,   // Bのヨコハバ
             LUp.devInfo);

  assert( status == CUSOLVER_STATUS_SUCCESS );

  double *X = (double*)malloc(sizeof(double)*LUp.n);
  // 結果を取得
  cudaMemcpy(X, LUp.dB, bytesof<double>(LUp.n*LUp.nrhs),
	     cudaMemcpyDeviceToHost);

  for(int i=0; i<LUp.n; i++) x[i] = X[i];


  return x;
}


static LUparameter LUinit(matrix<double>& AA)
{
  static LUparameter LUp;
  int i,j, n = AA.size();

  double *A = (double*)malloc(sizeof(double)*n*n);

  if (A == NULL) {
    fprintf(stderr,"can't malloc()\n");
    exit(0);
  }
  for(i=0;i<n*n;i++) A[i] = 0.0;


  for (long i=0; i<n; i++ ) {
    auto AAi = AA[i];
    auto j = AAi.begin();
    for ( ; j != AAi.end(); j++ )
      A[i*n+j->first] = AA[i][j->first];
  }
 
  int nrhs = 1;
  cusolverStatus_t status;
  
  // dense LAPACK
  printf("hello\n");
  status = cusolverDnCreate(&(LUp.handle));
  printf("world\n");
  assert( status == CUSOLVER_STATUS_SUCCESS );
  double* dA = allocate<double>(n*n);
  cudaMemcpy(dA, A, bytesof<double>(n*n), cudaMemcpyHostToDevice);
  int lda = n;

  // 必要なバッファ量を求め、確保する
  int worksize;
  status = cusolverDnDgetrf_bufferSize(
	     LUp.handle,
             n,   // 行
             n,   // 列
             dA,  // A
             lda, // Aのヨコハバ
             &worksize);
  assert( status == CUSOLVER_STATUS_SUCCESS );

  double* workspace = allocate<double>(worksize);

  // 計算結果に関する情報
  int* devInfo = allocate<int>(1);
  // ピボット
  int* pivot = allocate<int>(n);

  
  // LU分解 : dAに結果が求まる(それとpivot)
  status = cusolverDnDgetrf(
             LUp.handle,
             n,   // 行
             n,   // 列
             dA,  // A
             lda, // Aのヨコハバ
             workspace,
             pivot,
             devInfo);

  assert( status == CUSOLVER_STATUS_SUCCESS );

  double* dB = allocate<double>(n*nrhs);

  LUp.workspace = workspace;
  LUp.dA = dA;
  LUp.dB = dB;
  LUp.devInfo = devInfo;
  LUp.pivot = pivot;
  LUp.n = n;
  LUp.nrhs = nrhs;
  LUp.lda = lda;
  free(A);  A = NULL;
  return LUp;
}


static void LUdestroy(LUparameter LUp)
{
  cudaFree(LUp.workspace);       LUp.workspace = NULL;
  cudaFree(LUp.pivot);           LUp.pivot     = NULL;
  cudaFree(LUp.dA);              LUp.dA        = NULL;
  cudaFree(LUp.dB);              LUp.dB        = NULL;
  cudaFree(LUp.devInfo);         LUp.devInfo   = NULL;
  cusolverDnDestroy(LUp.handle); LUp.handle    = NULL;
  cudaDeviceReset();
}


int gpusolver(matrix<double> &A, vector<double> &x,const vector<double> &b)
{
  LUparameter LUp = LUinit(A);
  x = LUsolve(LUp,b);
  LUdestroy(LUp);
  return 0;
}
