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
/ OPENDX WRITERS UTILITIES                                                   /
/                                                                            /
/ #Version 0.1 Experimental:   30/08/2001 Alessandro Veneziani               /
/                                                                            /
/  I am sorry for my C++ at a very beginner level !!!!                       /
/                                                                            /
/ #Purpose: Container for subroutines for OPENDX                             /
/                                                                            /
/                                                                            /
/---------------------------------------------------------------------------*/
#ifndef _OPENDX_WRTRS_
#define _OPENDX_WRTRS_

#define OPENDX_VERSION 4.0 // to be verified
#define VRANK 1 //Vector Rank
#define SRANK 0 //Scalar Rank
#define EPS_DX 1.e-30 // a "sort" of machine epsilon due to a DX bug.
#ifdef TETRAHEDRON
#define CELL "tetrahedra" 
// to be completed with the other elements
#else
#define CELL "tetrahedra"
#endif

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

namespace LifeV
{
/*! This subroutine writes the header part of a OPENDX file (which is constant
  in a time-dependent but fixed-mesh problem). For a general problem, in
  particular for instance for a P2 on LinearTetra Mesh: in such a case
  we need to compute the coordinates of all the supplementary nodes.
  For "iso" cases (whenever we have ALL the coordinates of the nodes, it is
  recommended using the other subroutine with the same name but with a
  different signature.
*/
template <typename TheMesh, typename TheDof, typename TheFem>
void wr_opendx_header( std::string fname, const TheMesh& mesh, const TheDof& dof,
                       TheFem& fem, std::string nam_fe )
{

    std::ofstream ofile( fname.c_str() );

    ASSERT( ofile, "Error: Output file cannot be open" ); //
    UInt nldpe = fem.refFE.nbDofPerEdge;
    UInt nldpv = fem.refFE.nbDofPerVertex;
    UInt nldpf = fem.refFE.nbDofPerFace;
    UInt nldpV = fem.refFE.nbDofPerVolume;
    UInt nlv = dof.numLocalVertices();
    UInt nle = dof.numLocalEdges();
    UInt nlf = dof.numLocalFaces();
    UInt nV = mesh.numVolumes();
    UInt ne = mesh.numEdges();
    UInt nv = mesh.numVertices();
    UInt nf = mesh.numFaces();
    UInt nldof = nle * nldpe + nlv * nldpv + nlf * nldpf + nldpV;
    //  _totalDof=nV*nldpV+ne*nldpe+nv*nldpv+nf*nldpf;
    UInt num_points = dof.numTotalDof(); // discuss this with Luca
    UInt num_points_supp = dof.numTotalDof() - nv; // discuss this with Luca
    std::vector<Real> supp_x( num_points_supp, 0.0 );
    std::vector<Real> supp_y( num_points_supp, 0.0 );
    std::vector<Real> supp_z( num_points_supp, 0.0 );
    Real x, y, z;
    char quotes = '"';

    // Coordinates writing

    ofile << std::endl; // added by Diego
    ofile << "object " << quotes << "pos" << quotes << std::endl;
    ofile << "class array ";
    ofile << "type float ";
    ofile << "rank " << VRANK << " shape " << NDIM;
    ofile << " items " << num_points;
    ofile << " data follows" << std::endl;

    UInt i, ie, index, j;
    UInt gcount;
    UInt lcount;

    // Vertex based Dof: the coordinates are available from the Pont List
    for ( i = 0;i < nv;++i )
        ofile << mesh.pointList[ i ].x() << " " << mesh.pointList[ i ].y() << " "
        << mesh.pointList[ i ].z() << std::endl;

    // Now I store the coordinates of the supplementary nodes in a temporary vector
    // Edge Based Dof
    gcount = 0;
    lcount = nlv - 1;
    if ( nldpe > 0 )
    {
        for ( ie = 1; ie <= nV; ++ie )
        {
            fem.updateJac( mesh.volumeList( ie ) );
            for ( i = 1; i <= nle; ++i )
            {
                fem.coorMap( x, y, z, fem.refFE.xi( i + lcount ), fem.refFE.eta( i + lcount ),
                             fem.refFE.zeta( i + lcount ) );
                index = mesh.localEdgeId( ie, i );
                supp_x[ index - 1 ] = x;
                supp_y[ index - 1 ] = y;
                supp_z[ index - 1 ] = z;
            }
        }
        gcount += ne;
        lcount += nle;
    }
    // Face  Based Dof
    if ( nldpf > 0 )
    {
        for ( ie = 1; ie <= nV; ++ie )
        {
            fem.updateJac( mesh.volumeList( ie ) );
            for ( i = 1; i <= nlf; ++i )
            {
                fem.coorMap( x, y, z, fem.refFE.xi( i + lcount ), fem.refFE.eta( i + lcount ),
                             fem.refFE.zeta( i + lcount ) );
                index = mesh.localFaceId( ie, i ) + gcount;
                supp_x[ index - 1 ] = x;
                supp_y[ index - 1 ] = y;
                supp_z[ index - 1 ] = z;

            }
        }
        gcount += nf;
        lcount += nlf;
    }

    // Volume  Based Dof
    if ( nldpV > 0 )
    {
        for ( ie = 1; ie <= nV; ++ie )
        {
            fem.updateJac( mesh.volumeList( ie ) );
            //Alain (11/07/02) just one dof on volume ! change if there is more...
            for ( i = 1; i <= 1; ++i )
            {
                fem.coorMap( x, y, z, fem.refFE.xi( i + lcount ), fem.refFE.eta( i + lcount ),
                             fem.refFE.zeta( i + lcount ) );
                index = ie + gcount;
                supp_x[ index - 1 ] = x;
                supp_y[ index - 1 ] = y;
                supp_z[ index - 1 ] = z;
            }
        }
    }

    for ( i = 0;i < num_points_supp;++i )
    {
        if ( fabs( supp_x[ i ] ) < EPS_DX )
            supp_x[ i ] = 0.; // needed due to dx bugs
        if ( fabs( supp_y[ i ] ) < EPS_DX )
            supp_y[ i ] = 0.; // needed due to dx bugs
        if ( fabs( supp_z[ i ] ) < EPS_DX )
            supp_z[ i ] = 0.; // needed due to dx bugs
        ofile << supp_x[ i ] << " " << supp_y[ i ] << " " << supp_z[ i ] << std::endl;
    }

    ofile << std::endl;

    // connectivity
    ofile << "object " << quotes << "con" << quotes << std::endl;
    ofile << "class array ";
    ofile << "type int ";
    // Alain: this test was not working in case of using mixed fe
    //#if defined(LINEAR_P1)
    // it is replaced by:
    if ( nam_fe == "P1" )
        ofile << "rank " << VRANK << " shape " << nldof;
    else if ( nam_fe == "P1bubble" )
        ofile << "rank " << VRANK << " shape " << nldof - 1;
    //#elif defined(LINEAR_P2)
    else if ( nam_fe == "P2" )
        ofile << "rank " << VRANK << " shape " << nldof - 6;
    //#endif
    else if ( nam_fe == "P2tilde" )
        ofile << "rank " << VRANK << " shape " << nldof - 7;
    ASSERT( ( nam_fe == "P1" || nam_fe == "P2" || nam_fe == "P2tilde" ||
              nam_fe == "P1bubble" ),
            "output header defined only for finite element P1, P2 or P2tilde" );

    //#if defined(LINEAR_P1)
    if ( nam_fe == "P1" || nam_fe == "P1bubble" )
        ofile << " items " << nV;
    //#elif defined(LINEAR_P2)
    else if ( nam_fe == "P2" || nam_fe == "P2tilde" )
        ofile << " items " << nV * 8;
    //#endif
    ofile << " data follows" << std::endl;

    //#if defined(LINEAR_P1)
    if ( nam_fe == "P1" )
    {
        for ( i = 0;i < nV;++i )
        {
            for ( j = 0;j < nldof;++j )
                ofile << dof.localToGlobal( i + 1, j + 1 ) - 1 << " "; //damned (C vs) Fortran

            ofile << std::endl;
        }
    }

    else if ( nam_fe == "P1bubble" )
    {
        for ( i = 0;i < nV;++i )
        {
            for ( j = 0;j < nldof - 1;++j )
                ofile << dof.localToGlobal( i + 1, j + 1 ) - 1 << " "; //damned (C vs) Fortran

            ofile << std::endl;
        }
    }



    //#elif defined(LINEAR_P2)
    else if ( nam_fe == "P2" || nam_fe == "P2tilde" )
    {
        /*
          for (i=0;i<nV;++i){
          for (j=0;j<nldof-6;++j)
          ofile << dof.localToGlobal(i+1,j+1)-1 << " ";//damned (C vs) Fortran
         
          ofile << std::endl;
          }
         
        */

        int * MyCon = new int[ nldof ];

        for ( i = 0;i < nV;++i )
        {
            for ( j = 0;j < nldof;++j )
            {
                MyCon[ j ] = dof.localToGlobal( i + 1, j + 1 ) - 1; //damned (C vs) Fortran
            }
            ofile << MyCon[ 4 ] << " ";
            ofile << MyCon[ 0 ] << " ";
            ofile << MyCon[ 7 ] << " ";
            ofile << MyCon[ 6 ] << " ";
            ofile << std::endl;

            ofile << MyCon[ 6 ] << " ";
            ofile << MyCon[ 2 ] << " ";
            ofile << MyCon[ 9 ] << " ";
            ofile << MyCon[ 5 ] << " ";
            ofile << std::endl;

            ofile << MyCon[ 4 ] << " ";
            ofile << MyCon[ 1 ] << " ";
            ofile << MyCon[ 8 ] << " ";
            ofile << MyCon[ 5 ] << " ";
            ofile << std::endl;

            ofile << MyCon[ 7 ] << " ";
            ofile << MyCon[ 8 ] << " ";
            ofile << MyCon[ 9 ] << " ";
            ofile << MyCon[ 3 ] << " ";
            ofile << std::endl;

            ofile << MyCon[ 4 ] << " ";
            ofile << MyCon[ 9 ] << " ";
            ofile << MyCon[ 7 ] << " ";
            ofile << MyCon[ 8 ] << " ";
            ofile << std::endl;

            ofile << MyCon[ 4 ] << " ";
            ofile << MyCon[ 9 ] << " ";
            ofile << MyCon[ 7 ] << " ";
            ofile << MyCon[ 6 ] << " ";
            ofile << std::endl;

            ofile << MyCon[ 4 ] << " ";
            ofile << MyCon[ 9 ] << " ";
            ofile << MyCon[ 5 ] << " ";
            ofile << MyCon[ 8 ] << " ";
            ofile << std::endl;

            ofile << MyCon[ 4 ] << " ";
            ofile << MyCon[ 9 ] << " ";
            ofile << MyCon[ 5 ] << " ";
            ofile << MyCon[ 6 ] << " ";
            ofile << std::endl;

        } // end for (i = ...
        delete[] MyCon;

    }
    //#endif
    ofile << "attribute " << quotes << "element type" << quotes
    << " string " << quotes << CELL << quotes << std::endl;
    ofile << "attribute " << quotes << "ref" << quotes
    << " string " << quotes << "positions" << quotes << std::endl;

    ofile << std::endl;

}


/*! This subroutine considers the "iso" cases, whenever the points of the mesh
  actually coincide with the ones adopted for the computation. In this case,
  the coordinates of the nodes are ALL immediately available and there is no
  need of recomputing or rereading the other coordinates (which on the base of
  the McKoy algorithm are not explicitly computed so far (see fields...update))
*/
template <typename RegionMesh, typename Dof>
void wr_opendx_header( std::string fname, const RegionMesh& mesh, const Dof& dof )
{
    std::ofstream ofile( fname.c_str() );

    ASSERT( ofile, "Error: Output file cannot be open" ); //

    UInt nV = mesh.numVolumes();
    UInt nv = mesh.numVertices();
    UInt nldof = mesh.numLocalVertices();

    char quotes = '"';

    // Coordinates writing

    ofile << "object " << quotes << "pos" << quotes << std::endl;
    ofile << "class array ";
    ofile << "type float ";
    ofile << "rank " << VRANK << " shape " << NDIM;
    ofile << " items " << nv;
    ofile << " data follows" << std::endl;

    // Dof: the coordinates are available from the Pont List
    for (UInt i = 0;i < nv;++i )
    {
        ofile << mesh.pointList[ i ].x() << " " << mesh.pointList[ i ].y() << " "
        << mesh.pointList[ i ].z() << std::endl;
    }

    // connectivity
    ofile << "object " << quotes << "con" << quotes << std::endl;
    ofile << "class array ";
    ofile << "type int ";
    ofile << "rank " << VRANK << " shape " << nldof;
    ofile << " items " << nV;
    ofile << " data follows" << std::endl;

    for (UInt i = 0; i < nV; ++i )
    {
        for (UInt j = 0; j < nldof; ++j )
            ofile << dof.localToGlobal( i + 1, j + 1 ) - 1 << " "; //damned (C vs) Fortran

        ofile << std::endl;
    }

    ofile << "attribute " << quotes << "element type" << quotes << " string "
    << quotes << CELL << quotes << std::endl;
    ofile << "attribute " << quotes << "ref" << quotes << " string " << quotes
    << "positions" << quotes << std::endl;
    ofile << std::endl;


}

/*! {\tt  wr_opendx_scalar} is a subrotuine for writing SCALARS unknown for
  OPENDX. The file is appended to an existing one (at least you should call
  before {\tt wr_vtk_ascii_header}). Modification for accepting various type of
  vector. Alain 21/01/02.
*/
template <typename VectorType>
void wr_opendx_scalar( std::string fname, std::string name, VectorType const & U )
{
    unsigned int i;
    char quotes = '"';
    std::ofstream ofile( fname.c_str(), std::ios::app );
    std::string Ext = ".data"; // added by Diego
    std::string nameExt = name + Ext; // added by Diego

    ASSERT( ofile, "Error: Output file cannot be open" );


    ofile << "object " << quotes << nameExt << quotes << std::endl; // modif. by Diego
    ofile << "class array ";
    ofile << "type double ";
    ofile << "rank " << SRANK << " ";
    ofile << "items " << U.size() << " ";
    ofile << "data follows" << std::endl;

    for ( i = 0;i < U.size();++i )
        fabs(U[i]) < EPS_DX ? ofile << 0. << std::endl : ofile << U[i] << std::endl;

    ofile << "attribute " << quotes << "dep" << quotes << " string " << quotes
    << "positions" << quotes << std::endl;
    ofile << std::endl;
    ofile << "object " << quotes << name << quotes << " class field" << std::endl;
    ofile << "component " << quotes << "positions" << quotes << " value "
    << quotes << "pos" << quotes << std::endl;
    ofile << "component " << quotes << "connections" << quotes << " value "
    << quotes << "con" << quotes << std::endl;
    ofile << "component " << quotes << "data" << quotes << " value " << quotes
    << nameExt << quotes << std::endl; // modified by Diego
    ofile << std::endl; // added by Diego
    ofile << "end"; // added by Diego
    // caution: "end" must be put only once at the end of the file
    ofile << std::endl;
}

void wr_opendx_vector( std::string fname, std::string name, std::vector<Real> U,
                       std::vector<Real> V, std::vector<Real> W );

// Version accepting various vector types, Alain 21/01/02
template <typename VectorType>
void wr_opendx_vector( std::string fname, std::string name, VectorType U,
                       UInt nbcomp )
{
    unsigned int i, j, udim;
    char quotes = '"';
    std::ofstream ofile( fname.c_str(), std::ios::app );
    std::string Ext = ".data";
    std::string nameExt = name + Ext;

    ASSERT( ofile, "Error: Output file cannot be open" );

    udim = U.size();

    ofile << "object " << quotes << nameExt << quotes << std::endl;
    ofile << "class array ";
    ofile << "type float ";
    ofile << "rank " << VRANK << " shape " << nbcomp;
    ofile << " items " << udim / nbcomp << " ";
    ofile << "data follows" << std::endl;

    for ( i = 0;i < udim / nbcomp; ++i )
    {
        for ( j = 0; j < nbcomp; j++ )
        {
            unsigned int pos = i + j * udim / nbcomp;
            if ( fabs( U[ pos ] ) < EPS_DX )
                U[ pos ] = 0. ; // needed due to dx bugs
            ofile << U[ pos ] << " ";
        }
        ofile << std::endl;
    }

    ofile << "attribute " << quotes << "dep" << quotes << " string " << quotes
    << "positions" << quotes << std::endl;
    ofile << std::endl;
    ofile << "object " << quotes << name << quotes << " class field" << std::endl;
    ofile << "component " << quotes << "positions" << quotes << " value "
    << quotes << "pos" << quotes << std::endl;
    ofile << "component " << quotes << "connections" << quotes << " value "
    << quotes << "con" << quotes << std::endl;
    ofile << "component " << quotes << "data" << quotes << " value " << quotes
    << nameExt << quotes << std::endl;
    ofile << std::endl;
    ofile << "end";
    ofile << std::endl;
}
}
#endif
