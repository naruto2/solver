#ifndef OPERATOR_H
#define OPERATOR_H

#define M_solve(p) p
#define M_trans_solve(x) x
static int solveGMRES=0;

template < class Vector >
void printv(Vector &x)
{
  long i, n = x.size();
  for ( i = 0; i < n; i++ ) {
    cout<< x[i] << endl;
  }
}

void M_init(matrix<double>& A)
{
}

void M_destory(matrix<double>& A)
{
}


vector<double>& M_solveGMRES(vector<double>& x)
{
  return x;
}

vector<double>& operator+(vector<double>& y, vector<double>& x){
  int i, n = x.size();
  static vector<double> z(n);
  for ( i = 0; i < n; i++ ) z[i] = y[i] + x[i];
  return z;
}



vector<double>& operator-(vector<double>& y,  vector<double>& x) {
  int i, n = x.size();
  static vector<double> z(n);
  for( i = 0; i < n; i++ ) z[i] = y[i] - x[i];
  return z;
}


double dot(const vector<double>&x, const vector<double>&y){
  double rho = 0.0;
  long i, n = x.size();
  for ( i = 0; i < n; i++ ) rho += x[i]*y[i];
  return rho;
}

double nrm2(const vector<double>& r) {
  double norm = 0.0;
  long i, n = r.size();
  for ( i = 0; i < n; i++ ) norm += r[i]*r[i];
  return sqrt(norm);
}


vector<double>& operator*(const matrix<double>& A, vector<double>& b) {
  long n = b.size();
  static vector<double> x(n);
  for (long i=0; i<n; i++ ) {
    auto Ai = A[i];
    auto j = Ai.begin();
    for ( x[i]=0.0; j != Ai.end(); j++ )
      x[i] += j->second * b[j->first];
  }
  return x;
}
vector<double>& operator*(double a, vector<double>&x){
  long i, n = x.size();
  static vector<double> y(n);
  for ( i = 0; i < n; i++ ) y[i] = a*x[i];
  return y;
}


void cp(vector<double> &x, vector<double> &y){
  long i, n = x.size();
  for ( i = 0; i < n; i++ ) y[i] = x[i];
}

void axpy(long n, double *palpha,
	  const vector<double>& x, int l, vector<double>& y, int I) {
  long i;
  for ( i=0; i<n; i++) y[i] += *palpha*x[i];
}


template <class Vector>
Vector& y_ax(Vector& y, double alpha, const Vector& x) {
  axpy(x.size(), &alpha, x, 1, y, 1);
  return y;
}

vector<double>& trans_mult(matrix<double>& A, vector<double>& x)
{
  return A*x; //いずれ直す
}

template <class Matrix, class Vector>
void getProb(Matrix& A, Vector& x, Vector& b)
{
  long i, n = 1024;
  A.resize(n); x.resize(n); b.resize(n);
        // 中央差分行列を準備する
        //（対角項が2でその隣が1になる、↓こんなやつ）
        // | 2 1 0 0 0 0 0 0 ・・・ 0 0 0|
        // | 1 2 1 0 0 0 0 0 ・・・ 0 0 0|
        // | 0 1 2 1 0 0 0 0 ・・・ 0 0 0|
        // | 0 0 1 2 1 0 0 0 ・・・ 0 0 0|
        // | 0 0 0 1 2 1 0 0 ・・・ 0 0 0|
        // | 0 0 0 0 1 2 1 0 ・・・ 0 0 0|
        // | 0 0 0 0 0 1 2 1 ・・・ 0 0 0|
        // | 0 0 0 0 0 0 1 2 ・・・ 0 0 0|
        // | 0 0 0 0 0 0 0 0 ・・・ 2 1 0|
        // | 0 0 0 0 0 0 0 0 ・・・ 1 2 1|
        // | 0 0 0 0 0 0 0 0 ・・・ 0 1 2|

  for (i = 0; i < n; i++) {
    A[i][i] = 2;
    if(i > 0) A[i][i-1] = 1;
    if(i < n-1) A[i][i+1] = 1;
  }
  CRSinit(A);
  Vector x0(n), t;
  for (i = 0; i < n; i++) x0[i] = i*0.1;
  t =  A * x0;
  cp(t,b);
  CRSdestory(A);
}
#endif
