/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
#ifndef _USERFUNCTION_H
#define _USERFUNCTION_H
#include "user_diffusion.hpp"

namespace LifeV
{
double zero(const double& t, const double& x, const double& y, const double& z, const ID& i);
double g1(const double& t, const double& x, const double& y, const double& z, const ID& i);
double g2(const double& t, const double& x, const double& y, const double& z, const ID& i);
double g3(const double& t, const double& x, const double& y, const double& z, const ID& i);
double mixte_coeff(const double& t, const double& x, const double& y, const double& z, const ID& i);
double alpha(const double& t, const double& x, const double& y, const double& z, const ID& i);
double beta(const double& t, const double& x, const double& y, const double& z, const ID& i);
double p_adv(const double& t, const double& x, const double& y, const double& z, const ID& i);

class SourceFct
{
public:
  inline double operator()(double x,double y,double z,int ic=0) const {
    return 0.;
  }
};

//! Analytical solution (0 bc) and its Derivatives functions :
double Vfct(const double& t, const double& x, const double& y, const double& z, const ID& i); //!< Analytical solution

double VfctDer1(const double& t, const double& x, const double& y, const double& z, const ID& i); //! 1rst derivatives
double VfctDer2(const double& t, const double& x, const double& y, const double& z, const ID& i); //! 2nd derivatives
double minusLaplaceVfct(const double& t, const double& x, const double& y, const double& z, const ID& i); //!< := -laplace(Vfct)

class SourceAnalyticalFct
{
  // this is a check function: with this rhs, the solution
  // to the laplace pb should be (Ufct)
  // (with non homegeneous dirichlet b.c.)
public:
  inline double operator()(double x,double y,double z,int ic=0) const
  {
    return - minusLaplaceVfct(0., x, y, z, 1);
  }
};


class AnalyticalSolPres
{
public:
  inline double operator()(double x,double y,double z) const {
    //    return x*(1-x)*y*(1-y)*z*(1-z);
    return Vfct(0., x, y, z, 1);
  }
  inline double der_x(double x,double y,double z) const {
    //    return (1-2*x)*y*(1-y)*z*(1-z);
    return VfctDer1(0., x, y, z, 1);
  }
  inline double der_y(double x,double y,double z) const {
    //    return  x*(1-x)*(1-2*y)*z*(1-z);
    return VfctDer1(0., x, y, z, 2);
  }
  inline double der_z(double x,double y,double z) const {
    //    return x*(1-x)*y*(1-y)*(1-2*z);
    return VfctDer1(0., x, y, z, 3);
  }
};

class AnalyticalSolFlux
{
public:
  inline double operator()(double t,double x,double y,double z,int ic) const {
    switch(ic){
    case 1:
      return -( VfctDer1(0., x, y, z, 1) );
      break;
    case 2:
      return -( VfctDer1(0., x, y, z, 2) );
      break;
    case 3:
      return -( VfctDer1(0., x, y, z, 3) );
      break;
    default:
      std::cerr <<" AnalyticalSolFlux icomp ?????" << std::endl;
      exit(1);
    }
  }
};
}

#endif
