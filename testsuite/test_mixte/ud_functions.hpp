Real g(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return x*x + y*y + z*z;
    break;
  }
  return 0;
}

Real h(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return 2*z+2*(x*x + y*y + z*z);
  }
  return 0;
}


Real coef(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return 2;
    break;
  }
  return 0;
}
