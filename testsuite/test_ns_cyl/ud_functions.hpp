Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    return 0.;
}

Real u1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.0;
}

Real u3(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.1;
}

// Initial velocity 
Real u0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.0;
}

// Initial pressure
Real p0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
 return 1.0;
}

Real u2(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  Real pi = 3.14159265358979;  
  switch(i) {
  case 1:
  case 2:
    return 0.0;
    break;
  case 3:
     return 2.0*(sin(pi*t/2.0))*(1.0-4.0*(x*x+y*y));
     break; 
  }  
  return 0.0;
}
