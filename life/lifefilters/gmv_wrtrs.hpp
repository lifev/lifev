
/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
 
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
#ifndef _GMV_WRTRS_H
#define _GMV_WRTRS_H
#include <fstream>
#include "lifeV.hpp"

namespace LifeV
{
//using namespace std;
/*!
  \file gmv_wrtrs.hpp
  \author M.A. Fernandez
  \date 01/01/2004
 
  \brief Write a mesh and the Navier-Stokes solution for GMV Visualization Software
         http://www.mathematik.uni-dortmund.de/mirror/gmv/
 
 
*/

template <typename Mesh>
void wr_gmv_ascii( std::string fname, const Mesh& mesh, const UInt& dim, const Real* U, const Real* P )
{

    std::ofstream ofile( fname.c_str() );

    ASSERT( ofile, "Error: Output file cannot be open" );


    ofile << "gmvinput ascii\n";

    //
    // Nodes
    //
    ofile << "nodev ";
    UInt nV = mesh.numVertices();
    ofile << nV << std::endl;

    for ( UInt i = 1; i <= nV; ++i )
    {
        ofile << mesh.point( i ).x() << " "
        << mesh.point( i ).y() << " "
        << mesh.point( i ).z() << " "
        << std::endl;
    }
    ofile << std::endl;

    //
    // Elements
    //
    ofile << "cells ";
    UInt nE = mesh.numVolumes();
    ofile << nE << std::endl;

    typedef typename Mesh::VolumeShape ElementShape;
    UInt nVpE = ElementShape::numVertices;

    for ( ID k = 1; k <= nE; ++k )
    {
        ofile << "tet 4\n";
        for ( ID i = 1; i <= nVpE; ++i )
            ofile << mesh.volume( k ).point( i ).id() << " ";
        ofile << std::endl;
    }
    ofile << std::endl;
    //
    // Velocity
    //
    ofile << "velocity 1\n";
    for ( int icomp = 0; icomp < 3 ; ++icomp )
        for ( int i = 0; i < int( nV ); i++ )
        {
            ofile << U[ icomp * dim + i ] << std::endl;
        }

    ofile << std::endl;
    //
    // Pressure
    //
    ofile << "variable 'pressure' 1\n";
    for ( int i = 0; i < int( nV ); i++ )
    {
        ofile << P[ i ] << std::endl;
    }
    ofile << "endvars\n";
    ofile << std::endl;

    ofile << "endgmv\n";

}
}

#endif
