#ifndef USRFCT_H
#define USRFCT_H

#include "lifeV.hpp"

Real g1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return 1;
    break;
  }
  return 0;
}

Real g2(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return 2*z; //2*z;
    break;
  }
  return 0;
}

Real g3(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return -1; 
    break;
  }
  return 0;
}
#endif
