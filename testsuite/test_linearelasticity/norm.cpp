#include "norm.hpp"



Real maxnorm(const Vector& v) {
  Real max = abs(v[0]);
  for (UInt i=1; i<v.size(); ++i) {
    if ( abs(v[i]) > max )
      max =  abs(v[i]);
  }
  return max;
}


