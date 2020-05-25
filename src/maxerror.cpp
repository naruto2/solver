#include <cmath>
#include <iostream>
using namespace std;
#include <vector>
#include <map>
#include "matrix.hpp"


double maxerror(vector<double> r)
{
  long i;
  double max = 0.0;
  for(i=0;i<r.size();i++)if ( fabs(r[i])>fabs(max) ) max= r[i];
  return fabs(max);
}
