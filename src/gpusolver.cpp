#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cassert>
#include <vector>
#include <map>
#include "matrix.hpp"
using namespace std;


template<typename T>
inline size_t bytesof(unsigned int s) { return s * sizeof(T); }


template<typename T>
T* allocate(unsigned int size) {
  T* result = nullptr;
  cudaError_t status = cudaMalloc(&result, bytesof<T>(size));
  assert( status == cudaSuccess );
  return result;
}


static void Host2Device(double *h, double *d, unsigned int n)
{
  cudaMemcpy(h, d, bytesof<double>(n), cudaMemcpyHostToDevice);
}


static void Device2Host(double *h, double *d, unsigned int n)
{
  cudaMemcpy(h, d, bytesof<double>(n), cudaMemcpyDeviceToHost);
}


int gpusolver(matrix<double>&A, vector<double>&x, const vector<double> &b)
{
  static cusolverDnHandle_t handle;
  static cusolverStatus_t status;
  static double *a, *dA, *dB, *workspace;
  static int n, nrhs, lda, worksize, *pivot, *devInfo, init=1;

  n = A.size();
  nrhs    = 1;
  lda     = n; 
  dA      = allocate<double>(n*n);
  dB      = allocate<double>(n*nrhs);
  pivot   = allocate<int>(n);
  devInfo = allocate<int>(1);
  a       = (double*)calloc(n*n,sizeof(double));
  
  for (int i=0; i<n; i++ ) {
    auto Ai = A[i];
    auto j = Ai.begin();
    for ( ; j != Ai.end(); j++ )
      a[i*n+j->first] = A[i][j->first];
  }

  if (init) {
    status = cusolverDnCreate(&handle);
    printf("handle = %d init=%d\n",handle,init);
    init =0;
  }
  Host2Device(dA,a,n*n);

  status = cusolverDnDgetrf_bufferSize( handle, n, n, dA, lda, &worksize);
  printf("hello worksize =%d n=%d dA=%p lda=%d\n",worksize,n,dA,lda);
  workspace = allocate<double>(worksize);
  
  status = cusolverDnDgetrf( handle, n, n, dA, lda, workspace, pivot, devInfo);
  printf("world\n");
  Host2Device(dB, (double*)&b[0], n*nrhs);
  
  status = cusolverDnDgetrs( handle, CUBLAS_OP_T, n, nrhs, dA, lda, pivot,
			     dB, n, devInfo);
  printf("end\n");
  x.resize(n); 
  Device2Host((double*)&x[0], dB, n*nrhs);
  
  cudaFree(workspace);       
  cudaFree(dA);              
  cudaFree(dB);              
  cudaFree(devInfo);         
  cudaFree(pivot);
  //cusolverDnDestroy(handle);
  //cudaDeviceReset();
  free(a);
  return 0;
}
