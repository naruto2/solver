#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
#include <vector>
#include <map>
#include "matrix.hpp"
#include "operator.hpp"
#include "mmio.h"
#include "matrixmarket.h"




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
