#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
#include "matrix.hpp"
#include "operator.h"
#include "maxerror.h"
#include "mmio.h"
#include "matrixmarket.h"



void LU(matrix<double> &A, vector<double> &x, const vector<double> &b);


int main(int argc, char **argv)
{
  printf("%s ",argv[1]); fflush(stdout);

  matrix<double> A; vector<double> r, x, b;
  MatrixMarket(argv[1],A,x,b);

  LU(A,x,b);
  r = b - A*x;
  printf("%e\n",maxerror(r));
  return 0;
}
