Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    return 0.;
}

Real u1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.0;
}

Real fZero(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.0;
}

// Initial velocity 
Real u0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.0;
}


Real u2(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  //  Real pi = 3.14159265358979;  
  switch(i) {
  case 1:
  case 2:
    return 0.0;
    break;
  case 3:
    if (t <= 0.003 ) 
      //return 6500*( 1.0 - cos( pi*t / 2.5e-3 )  );
      return 1.3332e4;
    else 
      return 0.0;
    break; 
  }  
  return 0;
}


// Initial displacement and velocity 
Real d0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  switch(i) {
  case 1:
    //    return -z*(z-5)*x/50;
    return 0.;
    break;
  case 2:
    //    return -z*(z-5)*y/50;
    return 0.;
    break;
  case 3:
    return 0.;
    break;  
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break;
  }  
}

Real w0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
 
  switch(i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return 0.0;
    break;   
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break;
  }  
}


