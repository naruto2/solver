#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_
#include <map>
#include <vector>
using namespace std;

template<typename Real>
class matrix : public std::vector< std::map<long,Real> > {

public:
  matrix()       : std::vector< std::map<long, Real> >(){}
  matrix(long n) : std::vector< std::map<long, Real> >(n){}
};
#endif
