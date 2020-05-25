#ifndef MATRIXMARKET_H
#define MATRIXMARKET_H

template < class Matrix, class Vector >
  void MatrixMarket(char *argv1, Matrix& A, Vector& x, Vector &b)
{
  static double *val; static int *I, *J;
  int M, N, nz, ret;

  ret = mm_read_unsymmetric_sparse(argv1,&M,&N,&nz,&val,&I,&J);

  if ( ret != 0 ) printf("mm_read_unsymmetric_sparse()=%d %s\n",ret,
			 argv1);
  if ( M != N )return;

  A.resize(N); x.resize(N); b.resize(N);

  for (int k=0;k<nz;k++) A[I[k]][J[k]] = val[k];

  for (int k=0;k<N;k++) x[k] = 1.0;
  //CRSinit(A);
  {
    Vector t;
    b = A*x;
  }
  //CRSdestory(A);
  for (int k=0;k<N;k++) x[k] = 0.0;
}

#endif
