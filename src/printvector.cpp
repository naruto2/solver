#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <vector>
#include <map>
#include "matrix.hpp"
using namespace std;

void printvector(const vector<double>&b) {
  FILE *pp;

  pp = fopen("bar.txt","w");

  for (int i=0; i<b.size(); i++)
    fprintf(pp,"%e\n",b[i]);
  fflush(pp);
  getchar();
  fclose(pp);
}
