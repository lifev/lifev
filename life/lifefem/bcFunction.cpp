/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): M.A. Fernandez
             Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-12

  Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file bcFunction.cpp
   \author M.A. Fernandez
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-12
 */
#include <life.hpp>

#include <bcFunction.hpp>

namespace LifeV
{
//
// BCFunctionBase
//
BCFunctionBase::BCFunctionBase( const BCFunctionBase& bcf )
    :
    _M_g( bcf._M_g )
{
}

//! Constructor for BCFuncion_Base from a user defined function
BCFunctionBase::BCFunctionBase( function_type g )
    :
    _M_g( g )
{
}

//! set the function after having built it.
void
BCFunctionBase::setFunction( function_type g )
{
    _M_g = g;
}

//! Overloading function operator by calling attribut _g
Real
BCFunctionBase::operator() ( const Real& t, const Real& x, const Real& y,
                             const Real& z, const ID& icomp ) const
{
    return _M_g( t, x, y, z, icomp );
}

BCFunctionBase*
createBCFunctionBase( BCFunctionBase const* __bc )
{
    return new BCFunctionBase( ( BCFunctionBase const& )*__bc );
}
// register BCFunctionBase in factory for cloning
const bool __bcbase = FactoryCloneBCFunction::instance().registerProduct( typeid(BCFunctionBase), &createBCFunctionBase );

//
// BCFunctionMixte
//
BCFunctionMixte::BCFunctionMixte( const BCFunctionMixte& bcf )
    :
    BCFunctionBase( bcf ),
    _M_coef( bcf._M_coef )
{
}

BCFunctionMixte::BCFunctionMixte( function_type g, function_type coef )
    :
    BCFunctionBase( g ),
    _M_coef( coef )

{
}

//! set the functions after having built it.
void
BCFunctionMixte::setFunctions_Mixte( function_type g, function_type coef )
{
    setFunction( g );
    _M_coef = coef;
}

Real
BCFunctionMixte::coef( const Real& t, const Real& x, const Real& y,
                       const Real& z, const ID& icomp ) const
{
    return _M_coef( t, x, y, z, icomp );
}

BCFunctionBase*
createBCFunctionMixte( BCFunctionBase const* __bc )
{
    return new BCFunctionMixte( ( BCFunctionMixte const& )*__bc );
}
// register BCFunctionMixte in factory for cloning
const bool __bcmixte = FactoryCloneBCFunction::instance().registerProduct( typeid(BCFunctionMixte), &createBCFunctionMixte );




}






namespace LifeV
{




BCFunctionUDepBase::BCFunctionUDepBase(function_type g ):_M_g(g)
{
}
BCFunctionUDepBase::BCFunctionUDepBase(const BCFunctionUDepBase& bcf ):
				_M_g(bcf._M_g)
{
}
void
BCFunctionUDepBase::setFunction(function_type g)
{
  _M_g=g;
}
Real
BCFunctionUDepBase::operator()(const Real& t, const Real& x, const Real& y,
                        const Real& z, const ID& i, const Real& u ) const
{
  return _M_g(t,x,y,z,i,u);
}

BCFunctionUDepBase*
createBCFunctionUDep( BCFunctionUDepBase const* __bc )
{
  return new BCFunctionUDepBase( ( BCFunctionUDepBase const& )*__bc );
}
const bool __bcUDepBase = FactoryCloneBCFunctionUDep::instance().registerProduct(
		typeid(BCFunctionUDepBase), &createBCFunctionUDep );





BCFunctionUDepMixte::BCFunctionUDepMixte(function_type g,function_type coef):
				BCFunctionUDepBase(g),_M_coef(coef)
{
}

BCFunctionUDepMixte::BCFunctionUDepMixte(const BCFunctionUDepMixte& bcf):
				BCFunctionUDepBase(bcf),_M_coef(bcf._M_coef)
{
}
void
BCFunctionUDepMixte::setFunctions_Mixte(function_type g, function_type coef )
{
  setFunction(g);
  _M_coef=coef;
}
Real
BCFunctionUDepMixte::coef(const Real& t, const Real& x, const Real& y,
            const Real& z, const ID& icomp, const Real& u ) const
{
  return _M_coef( t, x, y, z, icomp, u );
}

BCFunctionUDepBase*
createBCFunctionUDepMixte( BCFunctionUDepBase const* __bc )
{
  return new BCFunctionUDepMixte( ( BCFunctionUDepMixte const& )*__bc );
}
const bool __bcUDepMixte = FactoryCloneBCFunctionUDep::instance().registerProduct(
		typeid(BCFunctionUDepMixte), &createBCFunctionUDepMixte );




}
