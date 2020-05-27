#include <cmath>
#include <vector>
#include <map>
#include "matrix.hpp"
#include "mmio.hpp"
using namespace std;


vector<double>& operator-(vector<double>& y,  vector<double>& x) {
  static vector<double> z;
  int i, n = x.size();
  z.resize(1);
  z.clear();
  z.resize(n);

  for( i = 0; i < n; i++ ) z[i] = y[i] - x[i];
  return z;
}


vector<double>& operator*(matrix<double>& A, vector<double>& b) {
  static vector<double> x;
  long n = A.size();
  x.resize(1);
  x.clear();
  x.resize(n);

  for (int i=0; i<n; i++ ) {
    auto Ai = A[i];
    auto j = Ai.begin();
    for ( x[i]=0.0; j != Ai.end(); j++ )
      x[i] += j->second * b[j->first];
  }
  return x;
}


static double maxerror(vector<double>&r)
{
  long i;
  double max = 0.0;
  for(i=0;i<r.size();i++)if ( fabs(r[i])>fabs(max) ) max= r[i];
  return fabs(max);
}


int main(int argc, char **argv)
{
  matrix<double> A; vector<double> x, b, r;
    
  for (int i=1; i<argc; i++) {
    printf("%s ",argv[i]); fflush(stdout);

    MatrixMarket(argv[i],A,x,b);
    printf("n = %d ",A.size()); fflush(stdout);

    solver(A,x,b);
    r = b - A*x;
    printf("%e\n",maxerror(r));
  }
  return 0;
}
