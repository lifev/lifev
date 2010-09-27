/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): A. Fumagalli  <alessio.fumagalli@mail.polimi.it>
            M. Kern       <michel.kern@inria.fr>
      Date: 2010-09

 Copyright (C) 2001-2006 EPFL, Politecnico di Milano, INRIA
 Copyright (C) 2006-2010 Politecnico di Milano

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

#ifndef _numericalFluxes_C_
#define _numericalFluxes_C_

#include <life/lifefem/NumericalFluxes.hpp>
#include <life/lifealg/brent.hpp>

namespace LifeV
{

// Constructor of the abstract class
AbstractNumericalFlux::
AbstractNumericalFlux ( const vectorFunction& physicalFlux,
                        const vectorFunction& firstDerivativePhysicalFlux ):
    M_physicalFlux                ( physicalFlux ),
    M_firstDerivativePhysicalFlux ( firstDerivativePhysicalFlux )
{

    CONSTRUCTOR( "AbstractNumericalFlux" );

} // constructor

// Destructor of the abstract class
AbstractNumericalFlux::
~AbstractNumericalFlux ()
{

    DESTRUCTOR( "AbstractNumericalFlux" );

} // destructor

// Create the fluxDotNormal function
boost::function<Real ( const Real& )>
AbstractNumericalFlux::
computeFunctionDotNormal ( const vectorFunction& function, const KN<Real>& normal, const Real& t,
                           const Real& x, const Real& y, const Real& z, const Real& plusMinus )
{

    scalarFunction functionDotNormalBound;

    // Bind the anonymousnamespace::fluxDotNormal function with known quantities.
    functionDotNormalBound = boost::bind( &functionDotNormal, _1, function,
                                          normal, t, x, y, z, plusMinus );

    return functionDotNormalBound;

} // computeFluxDotNormal

// ######################################################################### //

GodunovNumericalFlux::
GodunovNumericalFlux ( const vectorFunction& physicalFlux,
                       const vectorFunction& firstDerivativePhysicalFlux,
                       const Real& brentToll,
                       const UInt& brentMaxIter ):
    AbstractNumericalFlux::AbstractNumericalFlux  ( physicalFlux, firstDerivativePhysicalFlux ),
    M_brentToll                                   ( brentToll ),
    M_brentMaxIter                                ( brentMaxIter )
{

    CONSTRUCTOR( "GodunovNumericalFlux" );

} // constructor


Real
GodunovNumericalFlux::
operator() ( const Real& leftState, const Real& rightState, const KN<Real>& normal,
             const Real& t, const Real& x, const Real& y, const Real& z )
{

    // The value of the flux
    Real fluxValue( static_cast<Real>(0) );

    // The argmin or argmax of flux dot normal
    Real minMax( static_cast<Real>(0) );

    // The normal flux function
    scalarFunction normalFlux;

    if ( rightState > leftState )
    {
        // Create the function f \cdot n
        normalFlux = this->computeFunctionDotNormal( M_physicalFlux, normal, t, x, y, z, +1 );

        // Compute the argmin f \cdot n
        minMax = brent( normalFlux, leftState, rightState, M_brentToll, M_brentMaxIter );

        // Compute the flux value
        fluxValue = normalFlux( minMax );
    }
    else
    {
        // Create the function - f \cdot n
        normalFlux = this->computeFunctionDotNormal( M_physicalFlux, normal, t, x, y, z, -1 );

        // Compute the argmin - f \cdot n
        minMax = brent( normalFlux, leftState, rightState, M_brentToll, M_brentMaxIter );

        // Compute the flux value
        fluxValue = - normalFlux( minMax );
    }

    return fluxValue;

} // operator()

}

#endif //_numericalFluxes_C_
