#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <vector>
#include <map>
#include "matrix.hpp"


extern "C" {
  FILE *popen(const char *s, const char *mode);
  int pclose(FILE *f);
}


void printmatrix(matrix<double>&A) {
  FILE *pp;

  pp = fopen("foo.txt","w");
  fprintf(pp,"unset border\n");
  fprintf(pp,"unset xtics\n");
 
  fprintf(pp,"set xrange [-1:%d]\n",A.size());
  fprintf(pp,"set yrange [%d:-1]\n",A.size());


  for (int i=0; i<A.size(); i++ ) {
    auto Ai = A[i];
    auto j = Ai.begin();
    for ( ; j != Ai.end(); j++ )
      fprintf(pp,"set label \"%e\" at %d, %d;\n",A[j->first][i],j->first,i);
  }

  fprintf(pp,"plot '-' with lines title \"\"\n");
  fprintf(pp,"0 0\n");
  fprintf(pp,"%d %d\n",A.size()-1,A.size()-1);
  fprintf(pp,"e\n\n");
  fflush(pp);
  getchar();
  fclose(pp);
}
