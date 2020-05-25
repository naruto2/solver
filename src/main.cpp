#include <cmath>
#include <vector>
#include <map>
#include "matrix.hpp"
#include "mmio.hpp"
using namespace std;


vector<double>& operator-(vector<double>& y,  vector<double>& x) {
  int i, n = x.size();
  static vector<double> z(n);
  for( i = 0; i < n; i++ ) z[i] = y[i] - x[i];
  return z;
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


static double maxerror(vector<double> r)
{
  long i;
  double max = 0.0;
  for(i=0;i<r.size();i++)if ( fabs(r[i])>fabs(max) ) max= r[i];
  return fabs(max);
}


int main(int argc, char **argv)
{
  printf("%s ",argv[1]); fflush(stdout);

  matrix<double> A; vector<double> r, x, b;
  MatrixMarket(argv[1],A,x,b);
  solver(A,x,b);
  r = b - A*x;
  printf("%e\n",maxerror(r));
  return 0;
}
