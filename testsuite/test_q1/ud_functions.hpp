Real g(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return x*x + y*y + z*z;
    break;
  }
  return 0;
}
Real g1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return 0;
    break;
  }
  return 0;
}
Real g2(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return 2*z;
    break;
  }
  return 0;
}
Real g3(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return 10;
    break;
  }
  return 0;
}

Real h(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return 2*z;
  }
  return 0;
}


