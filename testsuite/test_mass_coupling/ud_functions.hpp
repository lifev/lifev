namespace LifeV
{
// **************** Function for fluid dynamics ********************************

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
  switch(i) {
  case 1:
  case 2:
    return 0.0;
    break;
  case 3:
     if(t <= 1.0)
	return 10.0*t*(1.0-4.0*(x*x+y*y));
     else	
	return 10.0*(1.0-4.0*(x*x+y*y));
     break; 
  }  
  return 0.0;
}

//  ************** Functions for mass transport ********************************

Real fc(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    return 0.;
}

Real c1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 1.0;
}

Real alpha(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 1.25513e-9;
}

Real beta(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
   return -1.75883e-6;
}

Real cw(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.3;
}

// Initial concentration in the lumen
Real cl_0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 1.0;
}

// Initial concentration in the wall
Real cw_0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.01;
}
// concentration at the adventitia
Real cadv(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.01;
}
}
