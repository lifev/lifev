// velocity imposed on the walls
Real u_in(const Real& t, const Real& x, const Real& y, const Real& z,
	const ID& i)
{
  switch(i){
  case 3:
     return 1.;
     break;
  case 1: case 2:
    return 0.;
    break;
  }
  return 0.;
}

//
Real u_w(const Real& t, const Real& x, const Real& y, const Real& z,
	const ID& i)
{
  return 0.;
}
//
Real u_out(const Real& t, const Real& x, const Real& y, const Real& z,
	const ID& i)
{
  return 0.;
}
