#ifndef SPARSE_H
#define SPARSE_H
#include <map>
#include <vector>

namespace sparse {
  using namespace std;

  template<typename Real> class matrix : public vector< map<long,Real> > {
  public:
    matrix()       : vector< map<long, Real> >(){}
    matrix(long n) : vector< map<long, Real> >(n){}
  };
}
#endif
