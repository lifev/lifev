/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
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
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-12
 */
#include <lifeV.hpp>

#include <bcFunction.hpp>

namespace LifeV
{
//
// BCFunctionBase
//
BCFunctionBase::BCFunctionBase( const BCFunctionBase& bcf )
{
    _g = bcf._g;
}

//! Constructor for BCFuncion_Base from a user defined function
BCFunctionBase::BCFunctionBase( Function g )
{
    _g = g;
}

//! set the function after having built it.
void BCFunctionBase::setFunction( Function g )
{
    _g = g;
}

//! Overloading function operator by calling attribut _g
Real
BCFunctionBase::operator() ( const Real& t, const Real& x, const Real& y,
                             const Real& z, const ID& icomp ) const
{
    return _g( t, x, y, z, icomp );
}

//
// BCFunctionMixte
//
BCFunctionMixte::BCFunctionMixte( const BCFunctionMixte& bcf )
    :
    BCFunctionBase( bcf._g )
{
    _coef = bcf._coef;
}

BCFunctionMixte::BCFunctionMixte( Function g, Function coef )
    :
    BCFunctionBase( g )
{
    _coef = coef;
}

//! set the functions after having built it.
void BCFunctionMixte::setFunctions_Mixte( Function g, Function coef )
{
    _g = g;
    _coef = coef;
}

Real BCFunctionMixte::coef( const Real& t, const Real& x, const Real& y,
                             const Real& z, const ID& icomp ) const
{
    return _coef( t, x, y, z, icomp );
}


}

