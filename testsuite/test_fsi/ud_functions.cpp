/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file ud_functions.cpp
*/

#include "ud_functions.hpp"


namespace LifeV
{
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
  switch(i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
     if ( t <= 0.003 )
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
    return 0.;
    break;
  case 2:
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
}

