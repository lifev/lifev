/*
 This file is part of the LifeV library
 Copyright (C) 2004 EPFL, INRIA and Politecnico di Milano

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
  \file dataConvDiffReact.hpp
  \author M. Prosi
  \date 03/2004
  \version 1.0

  \brief File containing a class for handling Convection-Diffusion-Reaktion processes data with GetPot

*/
#ifndef _DATACONVDIFFREACT_H_
#define _DATACONVDIFFREACT_H_
#include <string>
#include <iostream>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>

namespace LifeV
{
/*!
  \class DataConvDiffReact

  Base class which holds usual data for the Convection-Diffusion-Reaction equation solvers

*/
template <typename Mesh>
class DataConvDiffReact: public DataMesh<Mesh>, public DataTime
{
public:

    //! Constructor
    DataConvDiffReact( const GetPot& dfile );

    //! Ouptut
    void showMe( std::ostream& c = std::cout );
    //! Diffusivity
    Real diffusivity() const;
    //! Reaction coefficient
    Real react() const;
    int stationary() const;

protected:
    //! Physics
    Real _diffusivity; // Diffusivity
    Real _react; // Reaction coefficient
    int _stationary; // switch stationary/instationary calculation
};


//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh>
DataConvDiffReact<Mesh>::
DataConvDiffReact( const GetPot& dfile ) :
        DataMesh<Mesh>( dfile, "masstransport/discretization" ),
        DataTime( dfile, "masstransport/time" )
{

    // physics
    _diffusivity = dfile( "masstransport/physics/diffusivity", 1. );
    _react = dfile( "masstransport/physics/react", 1. );
    _stationary = dfile( "masstransport/physics/stationary", 0 );
}

// Output
template <typename Mesh>
void DataConvDiffReact<Mesh>::
showMe( std::ostream& c )
{
    // physics
    c << "\n*** Values for data [masstransport/physics]\n\n";
    c << "diffusivity   = " << _diffusivity << std::endl;
    c << "reaction coefficient  = " << _react << std::endl;
    c << "stationary = " << _stationary << std::endl;

    c << "\n*** Values for data [masstransport/discretization]\n\n";
    DataMesh<Mesh>::showMe( c );
    c << "\n*** Values for data [masstransport/time]\n\n";
    DataTime::showMe( c );

}
////////////////////
// The diffusivity
template <typename Mesh>
Real DataConvDiffReact<Mesh>::
diffusivity() const
{
    return _diffusivity;
}

// The reaction coefficient
template <typename Mesh>
Real DataConvDiffReact<Mesh>::
react() const
{
    return _react;
}

// The switch stationary/unstationary
template <typename Mesh>
int DataConvDiffReact<Mesh>::
stationary() const
{
    return _stationary;
}
}
#endif
