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
/*! --------------------------------------------------------------------------*
/                                                                            /
/      ...                                                                   /
/                                                                            /
/                                                                            /
/ GMSH WRITERS UTILITIES                                                     /
/                                                                            /
/ #Version 0.1 Experimental:   19/03/2001 Alessandro Veneziani               /
/                                                                            /
/  I am sorry for my C++ at a very beginner level !!!!                       /
/                                                                            /
/ #Purpose: Container for subroutines for GMSH                                /
/                                                                            /
/                                                                            /
/---------------------------------------------------------------------------*/
#ifndef _GMSH_WRTRS_
#define _GMSH_WRTRS_

#define GMSH_VERSION 1.0

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

namespace LifeV
{
template <typename TheMesh, typename TheDof, typename TheFem, typename GeoMap>
void wr_gmsh_parsed( std::string fname, std::string title, const TheMesh& mesh,
                     const TheDof& dof, const TheFem& fem, const GeoMap& geo,
                     std::vector<Real> U )
{
    std::ofstream ofile( fname.c_str() );

    ASSERT( ofile, "Error: Output file cannot be open" ); //

    UInt nldpe = fem.nbDofPerEdge;
    UInt nldpv = fem.nbDofPerVertex;
    UInt nldpf = fem.nbDofPerFace;
    UInt nldpV = fem.nbDofPerVolume;
    UInt nlv = dof.nlv();
    UInt nle = dof.nle();
    UInt nlf = dof.nlf();
    UInt nV = mesh.numVolumes();
    UInt ne = mesh.numEdges();
    UInt nv = mesh.numVertices();
    UInt nf = mesh.numFaces();
    UInt nldof = nle * nldpe + nlv * nldpv + nlf * nldpf + nldpV;
    //  _totalDof=nV*nldpV+ne*nldpe+nv*nldpv+nf*nldpf;
    UInt cells_size;
    UInt num_points = dof.numTotalDof(); // discuterne con Luca
    UInt num_points_supp = dof.numTotalDof() - nv; // discuterne con Luca
    std::vector<Real> supp_x( num_points_supp, 0.0 );
    std::vector<Real> supp_y( num_points_supp, 0.0 );
    std::vector<Real> supp_z( num_points_supp, 0.0 );
    Real x, y, z;
    char virgole = '"';

    ofile << "View " << virgole << title << virgole << "{" << std::endl;

    UInt i, ie, index, j;
    UInt gcount;
    UInt lcount;

    // Vertex based Dof: the coordinates are available from the Pont List
    for ( i = 0;i < nv;++i )
        ofile << "SP(" << mesh.pointList[ i ].x() << "," << mesh.pointList[ i ].y() << "," << mesh.pointList[ i ].z() << ")" << "{" << U[ i ] << "};" << std::endl;

    // Now I store the coordinates of the supplementary nodes in a temporary vector
    // Edge Based Dof
    gcount = 0;
    lcount = nlv - 1;
    if ( nldpe > 0 )
    {
        for ( ie = 1; ie <= nV; ++ie )
        {
            geo.update( mesh.volumeList( ie ), 1 );
            fem.update( mesh.volumeList( ie ), 1 );
            for ( i = 1; i <= nle; ++i )
            {
                geo.coorMap( x, y, z, fem.xi( i + lcount ), fem.eta( i + lcount ), fem.zeta( i + lcount ) );
                index = mesh.localEdgeId( ie, i );
                supp_x[ index - 1 ] = x;
                supp_y[ index - 1 ] = y;
                supp_z[ index - 1 ] = z;
            }
        }
    }
    // Face  Based Dof
    gcount += ne;
    lcount += nle;
    if ( nldpf > 0 )
    {
        for ( ie = 1; ie <= nV; ++ie )
        {
            geo.update( mesh.volumeList( ie ), 1 );
            fem.update( mesh.volumeList( ie ), 1 );
            for ( i = 1; i <= nlf; ++i )
            {
                geo.coorMap( x, y, z, fem.xi( i + lcount ), fem.eta( i + lcount ), fem.zeta( i + lcount ) );
                index = mesh.localFaceId( ie, i ) + gcount;
                supp_x[ index - 1 ] = x;
                supp_y[ index - 1 ] = y;
                supp_z[ index - 1 ] = z;

            }
        }
    }

    // Volume  Based Dof
    gcount += nf;
    lcount += nlf;
    if ( nldpV > 0 )
    {
        for ( ie = 1; ie <= nV; ++ie )
        {
            geo.update( mesh.volumeList( ie ), 1 );
            fem.update( mesh.volumeList( ie ), 1 );
            for ( i = 1; i <= nV; ++i )
            {
                geo.coorMap( x, y, z, fem.xi( i + lcount ), fem.eta( i + lcount ), fem.zeta( i + lcount ) );
                index = ie + gcount;
                supp_x[ index - 1 ] = x;
                supp_y[ index - 1 ] = y;
                supp_z[ index - 1 ] = z;
            }
        }
    }

    for ( i = nv;i < nv + num_points_supp;++i )
    {
        ofile << "SP(" << supp_x[ i - nv ] << "," << supp_y[ i - nv ] << "," << supp_z[ i - nv ] << ")" << "{" << U[ i ] << "};" << std::endl;
    }

    for ( ie = 1; ie <= nV; ++ie )
    {
        ofile << "SS(";

    }

    ofile << std::endl;
    ofile << "};" << std::endl;
}
}
#endif

