Real g1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return t*t*(x*x+y*y+z*z);//1.5;
    break;
  }
  return t*t*(x*x+y*y+z*z);
}

Real g2(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
  case 1:
    return t*t*(x*x+y*y+z*z);//5.1;
    break;
  }
  return t*t*(x*x+y*y+z*z);
}

Real nu(const Real& t){
  return 1./(t*t);
}

Real sigma(const Real& t){
  return -2./t;
}

Real u0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  switch(i){
    case 1:
     return t*t*(x*x+y*y+z*z);//5.1;
     break;
  }	
  return t*t*(x*x+y*y+z*z);
}

