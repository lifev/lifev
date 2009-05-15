/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
/*!
  \file dataSimplifiedStructure.h
  \author M.A. Fernandez
  \date 01/2003
  \version 1.0

  \brief File containing a class for handling data for reduced
         structural models (algegraic law and independent ring model)

*/

#ifndef _DATASIMPLIFIEDSTRUCTURE_H_
#define _DATASIMPLIFIEDSTRUCTURE_H_
#include <string>
#include <iostream>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>

namespace LifeV
{
/*!
  \class DataNavierStokes

  Base class which holds usual data for reduced
  structural models (algegraic law and independent ring model)

*/
template <typename Mesh>
class DataSimplifiedStructure:
            public DataMesh<Mesh>
{
public:

    //! Constructor
    DataSimplifiedStructure( const GetPot& dfile );

    //! Ouptut
    void showMe( std::ostream& c = std::cout ) const;

protected:
    //! Physics
    Real _rho; // densisty
    Real _h;  // thickness
    Real _E;  // Young modulus
    Real _nu; // Poisson coeficient
    Real _R0; // Radius
};



//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh>
DataSimplifiedStructure<Mesh>::
DataSimplifiedStructure( const GetPot& dfile ) :
        DataMesh<Mesh>( dfile, "solid/discretization" )
{
    // physics
    _rho = dfile( "solid/physics/density", 1. );
    _h = dfile( "solid/physics/thickness", 1. );
    _E = dfile( "solid/physics/young", 1. );
    _nu = dfile( "solid/physics/poisson", 0.5 );
    _R0 = dfile( "solid/physics/radius", 1. );
}


// Output
template <typename Mesh>
void DataSimplifiedStructure<Mesh>::
showMe( std::ostream& c ) const
{
    // physics
    c << "\n*** Values for data [solid/physics]\n\n";
    c << "density   = " << _rho << std::endl;
    c << "thickness = " << _h << std::endl;
    c << "young     = " << _E << std::endl;
    c << "poisson   = " << _nu << std::endl;
    c << "radius    = " << _R0 << std::endl;
    c << "\n*** Values for data [solid/discretization]\n\n";
    DataMesh<Mesh>::showMe();
}
}
#endif
