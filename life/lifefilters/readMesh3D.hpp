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
#ifndef _READMESH3D_HH_
#define _READMESH3D_HH_
#include "regionMesh3D.hpp"
#include "util_string.hpp"
#include "mesh_util.hpp"
#include <sstream>

namespace LifeV
{
// IT NEEDS THE NEW mesh_util.h (V 0.2 onwards)
// #define OLDMPPFILE // UNCOMMENT IF YOU HAVE AN OLD Mesh++ file!
/*----------------------------------------------------------------------*
|
|
|
|
|
| #Version 0.2 Experimental 19/8/99. Luca Formaggia
|
| Added support for numbering from 1
! Added markers
|
| #Mesh readers
|
*----------------------------------------------------------------------*/
/********************************************************************************
MESH BUILDERS
********************************************************************************/

bool
readMppFileHead( std::ifstream & mystream, UInt & numVertices, UInt & numBVertices, UInt & numBFaces, UInt & numBEdges, UInt & numVolumes );

//
//=================================================================
// ReadmppFile It reads mesh++ Tetra meshes. It convert them into
// Quadratic Tetra if needed it return false if  mesh check is uncessessfull
//
template <typename GeoShape, typename MC>
bool
readMppFile( RegionMesh3D<GeoShape, MC> & mesh,
             const std::string & filename,
             EntityFlag regionFlag,
             bool verbose=false )
{
    unsigned done = 0;
    std::string line;
    Real x, y, z;
    int ity, ity_id;
    UInt p1, p2, p3, p4;
    UInt nVe( 0 ), nBVe( 0 ), nFa( 0 ), nBFa( 0 ), nPo( 0 ), nBPo( 0 );
    UInt nVo( 0 ), nEd( 0 ), nBEd( 0 );
    UInt i;

    std::stringstream discardedLog;
    std::ostream& oStr = verbose ? std::cout : discardedLog;

    ASSERT_PRE0( GeoShape::Shape == TETRA , "ReadMppFiles reads only tetra meshes" ) ;

    ASSERT_PRE0( GeoShape::Shape == TETRA, "Sorry, ReadMppFiles reads only tetra meshes" );

    ASSERT_PRE0( GeoShape::numVertices <= 6, "Sorry, ReadMppFiles handles only liner&quad tetras" );

    // open stream to read header

    std::ifstream hstream( filename.c_str() );
    if ( hstream.fail() )
    {
        std::cerr << " Error in readMpp: File not found or locked" << std::endl;
        abort();
    }
    oStr << "Reading Mesh++ file" << std::endl;
    if ( ! readMppFileHead( hstream, nVe, nBVe, nBFa, nBEd, nVo ) )
    {
        std::cerr << " Error While reading mesh++ file headers" << std::endl;
        ABORT() ;
    }
    hstream.close();

    //Reopen the stream: I know it is stupid but this is how it goes
    std::ifstream mystream( filename.c_str() );
    if ( mystream.fail() )
    {
        std::cerr << " Error in readMpp: File not found or locked" << std::endl;
        abort();
    }

    // Euler formulas
    nFa = 2 * nVo + ( nBFa / 2 );
    nEd = nVo + nVe + ( 3 * nBFa - 2 * nBVe ) / 4;

    // Be a little verbose
    if ( GeoShape::numPoints > 4 )
    {

        std::cout << "Quadratic Tetra  Mesh (from Linear geometry)"
                  << std::endl;
        nPo = nVe + nEd;
        nBPo = nBVe + nBEd;
    }
    else
    {
        std::cout << "Linear Tetra Mesh" << std::endl;
        nPo = nVe;
        nBPo = nBVe;
    }
    std::cout << "#Vertices = "          << std::setw(10) << nVe
              << "  #BVertices       = " << std::setw(10) << nBVe << std::endl;
    std::cout << "#Faces    = "          << std::setw(10) << nFa
              << "  #Boundary Faces  = " << std::setw(10) << nBFa << std::endl;
    std::cout << "#Edges    = "          << std::setw(10) << nEd
              << "  #Boundary Edges  = " << std::setw(10) << nBEd << std::endl;
    std::cout << "#Points   = "          << std::setw(10) << nPo
              << "  #Boundary Points = " << std::setw(10) << nBPo << std::endl;
    std::cout << "#Volumes  = "          << std::setw(10) << nVo  << std::endl;

    // Set all basic data structure

    // I store all Points
    mesh.setMaxNumPoints( nPo, true );
    mesh.setNumBPoints( nBPo );
    mesh.numVertices() = nVe;
    mesh.numBVertices() = nBVe;
    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges( nBEd );
    mesh.numEdges() = nEd; // Here the REAL number of edges (all of them)
    mesh.setNumBEdges( nBEd );
    // Only Boundary Faces
    mesh.setMaxNumFaces( nBFa );
    mesh.numFaces() = nFa; // Here the REAL number of edges (all of them)
    mesh.setNumBFaces( nBFa );

    mesh.setMaxNumVolumes( nVo, true );

    mesh.setMarker( regionFlag ); // Mark the region

    typename RegionMesh3D<GeoShape, MC>::PointType * pp = 0;
    typename RegionMesh3D<GeoShape, MC>::EdgeType * pe = 0;
    typename RegionMesh3D<GeoShape, MC>::FaceType * pf = 0;
    typename RegionMesh3D<GeoShape, MC>::VolumeType * pv = 0;
    // addPoint()/Face()/Edge() returns a reference to the last stored point
    // I use that information to set all point info, by using a pointer.
    UInt count = 0;
    long int ibc;
    while ( next_good_line( mystream, line ).good() )
    {
        if ( line.find( "odes" ) != std::string::npos )
        {
            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );
            //      _numVertices=atoi(node_s);
            for ( i = 0;i < nVe;i++ )
            {
#ifdef OLDMPPFILE
                mystream >> x >> y >> z >> ity >> ibc;
#else

                mystream >> x >> y >> z >> ity >> ity_id;
                if ( ity != 3 )
                    mystream >> ibc;
#endif

                if ( ity != 3 )
                {
                    ++count;
                    pp = &mesh.addPoint( true ); // Boundary point. Boundary switch set by the mesh method.
                    pp->setMarker( EntityFlag( ibc ) );
                }
                else
                {
                    pp = &mesh.addPoint( false );
                }
                pp->x() = x;
                pp->y() = y;
                pp->z() = z;
                pp->setMarker( EntityFlag( ibc ) );
            }
            oStr << "Vertices Read " << std::endl;
            done++;
            if ( count != nBVe )
                std::cerr << "NumB points inconsistent !" << std::endl;
        }
        if ( line.find( "iangular" ) != std::string::npos )
        {
            oStr << "Reading Bfaces " << std::endl;
            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );
            // _numBFaces=atoi(node_s);
            for ( i = 0;i < nBFa;i++ )
            {
#ifdef OLDMPPFILE
                mystream >> p1 >> p2 >> p3 >> ity >> ibc;
#else

                mystream >> p1 >> p2 >> p3 >> ity >> ity_id >> ibc;
#endif

                pf = &( mesh.addFace( true ) ); // Only boundary faces

                pf->setMarker( EntityFlag( ibc ) );
                pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.
            }
            oStr << "Boundary Faces Read " << std::endl;
            done++;
        }
        if ( line.find( "Sides" ) != std::string::npos )
        {
            oStr << "Reading Bedges " << std::endl;
            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );
            //_numBEdges=atoi(node_s);
            for ( i = 0;i < nBEd;i++ )
            {
#ifdef OLDMPPFILE
                mystream >> p1 >> p2 >> ity >> ibc;
#else

                mystream >> p1 >> p2 >> ity >> ity_id >> ibc;
#endif

                pe = &mesh.addEdge( true ); // Only boundary edges.
                pe->setMarker( EntityFlag( ibc ) );
                pe->setPoint( 1, mesh.point( p1 ) ); // set edge conn.
                pe->setPoint( 2, mesh.point( p2 ) ); // set edge conn.
            }
            oStr << "Boundary Edges Read " << std::endl;
            done++;
        }
        count = 0;
        if ( line.find( "etrahedral" ) != std::string::npos )
        {
            oStr << "Reading Volumes " << std::endl;
            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );
            for ( i = 0; i < nVo; i++ )
            {
                mystream >> p1 >> p2 >> p3 >> p4;
                pv = &mesh.addVolume();
                pv->id() = i + 1;
                pv->setPoint( 1, mesh.point( p1 ) );
                pv->setPoint( 2, mesh.point( p2 ) );
                pv->setPoint( 3, mesh.point( p3 ) );
                pv->setPoint( 4, mesh.point( p4 ) );
                count++;
            }
            oStr << count << " Volume elements Read" << std::endl;
            done++;
        }

    }
    // This part is to build a P2 mesh from a P1 geometry

    if ( GeoShape::numPoints > 4 )
        p1top2( mesh );
    mystream.close();

    // Test mesh
    Switch sw;

    ///// CORRECTION JFG
    //if (mesh.check(1, true,true))done=0;

        if ( !checkMesh3D( mesh, sw, true, verbose, oStr, std::cerr, oStr ) )
        abort(); // CORRECTION JFG

    Real vols[ 3 ];
    getVolumeFromFaces( mesh, vols, oStr );
    oStr << "   VOLUME ENCLOSED BY THE MESH COMPUTED BY INTEGRATION ON" <<
    " BOUNDARY FACES" << std::endl;
    oStr << "INT(X)     INT(Y)      INT(Z) <- they should be equal and equal to" << std::endl <<
    "                                 the voulume enclosed by the mesh " << std::endl;
    oStr << vols[ 0 ] << " " << vols[ 1 ] << " " << vols[ 2 ] << std::endl;

    oStr << "   BOUNDARY FACES ARE DEFINING A CLOSED SURFACE IF " << testClosedDomain( mesh, oStr ) << std::endl <<
    " IS (ALMOST) ZERO" << std::endl;

    return done == 4 ;

}

/*
                                          INRIA MESH FILE READERS
       29/06/2002 L.F.
 */
//! INRIAMesh used either spaces or CR as separators
/*! It sucks, but this is the way it is! This little function should help handling it. It gets an
  integer field from the std::string line if it is not empty, otherwise from the input stream.

  It assumes that the std::string is either empty or it contains and integer!!! No check is made to verify this.
*/

int nextIntINRIAMeshField( std::string const & line, std::istream & mystream );

bool
readINRIAMeshFileHead( std::ifstream & mystream,
                       UInt & numVertices,
                       UInt & numBVertices,
                       UInt & numBFaces,
                       UInt & numBEdges,
                       UInt & numVolumes,
                       ReferenceShapes & shape );

//=================================================================
// ReadmppFile It reads mesh++ Tetra meshes. It convert them into
// Quadratic Tetra if needed it return false if  mesh check is uncessessfull
//
template <typename GeoShape, typename MC>
bool
readINRIAMeshFile( RegionMesh3D<GeoShape, MC> & mesh,
                   std::string const & filename,
                   EntityFlag regionFlag,
                   bool verbose=false )
{
    unsigned done = 0;
    std::string line;
    Real x, y, z;
    UInt p1, p2, p3, p4, p5, p6, p7, p8;
    UInt nVe( 0 ), nBVe( 0 ), nFa( 0 ), nBFa( 0 ), nPo( 0 ), nBPo( 0 );

    UInt nVo( 0 ), nEd( 0 ), nBEd( 0 );
    UInt i;
    ReferenceShapes shape;

    std::stringstream discardedLog;
    std::ostream& oStr = verbose ? std::cout : discardedLog;

    // open stream to read header

    std::ifstream hstream( filename.c_str() );
    if ( hstream.fail() )
    {
        std::cerr << " Error in readINRIAMeshFile: File not found or locked"
                  << std::endl;
        abort();
    }
    oStr << "Reading INRIA mesh file" << std::endl;
    if ( ! readINRIAMeshFileHead( hstream, nVe, nBVe, nBFa, nBEd, nVo, shape ) )
    {
        std::cerr << " Error While reading INRIA mesh file headers" << std::endl;
        ABORT() ;
    }
    hstream.close();

    //Reopen the stream: I know it is stupid but this is how it goes
    std::ifstream mystream( filename.c_str() );
    if ( mystream.fail() )
    {
        std::cerr << " Error in readINRIAMeshFile: File not found or locked" << std::endl;
        abort();
    }

    ASSERT_PRE0( GeoShape::Shape == shape, "INRIA Mesh file and mesh element shape is not consistent" );

    // Euler formulas to get number of faces and number of edges
    nFa = 2 * nVo + ( nBFa / 2 );
    nEd = nVo + nVe + ( 3 * nBFa - 2 * nBVe ) / 4;

    // Be a little verbose
    switch ( shape )
    {

    case HEXA:
        ASSERT_PRE0( GeoShape::numPoints == 8, "Sorry I can read only bilinear Hexa meshes" );
        oStr << "Linear Hexa Mesh" << std::endl;
        nPo = nVe;
        nBPo = nBVe;
        break;
    case TETRA:
        if ( GeoShape::numPoints > 4 )
        {
            //    if (GeoShape::numPoints ==6 ){
            std::cout << "Quadratic Tetra  Mesh (from Linear geometry)" << std::endl;
            nPo = nVe + nEd;
            // nBPo=nBVe+nBEd; // FALSE : nBEd is not known at this stage in a INRIA file (JFG 07/2002)
            // I use the relation  nBVe + nBFa - 2 = nBEd, But, is it general (hole...) ???? (JFG 07/2002)
            nBPo = nBVe + ( nBVe + nBFa - 2 );
        }
        else
        {
            std::cout << "Linear Tetra Mesh" << std::endl;
            nPo = nVe;
            nBPo = nBVe;
        }
        break;
    default:
        ERROR_MSG( "Current version of INRIA Mesh file reader only accepts TETRA and HEXA" );
    }

    std::cout << "#Vertices = "          << std::setw(10) << nVe
              << "  #BVertices       = " << std::setw(10) << nBVe << std::endl;
    std::cout << "#Faces    = "          << std::setw(10) << nFa
              << "  #Boundary Faces  = " << std::setw(10) << nBFa << std::endl;
    std::cout << "#Edges    = "          << std::setw(10) << nEd
              << "  #Boundary Edges  = " << std::setw(10) << nBEd << std::endl;
    std::cout << "#Points   = "          << std::setw(10) << nPo
              << "  #Boundary Points = " << std::setw(10) << nBPo << std::endl;
    std::cout << "#Volumes  = "          << std::setw(10) << nVo  << std::endl;

    // Set all basic data structure

    // I store all Points
    mesh.setMaxNumPoints( nPo, true );
    mesh.setNumBPoints( nBPo );
    mesh.numVertices() = nVe;
    mesh.numBVertices() = nBVe;
    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges( nBEd );
    mesh.numEdges() = nEd; // Here the REAL number of edges (all of them)
    mesh.setNumBEdges( nBEd );
    // Only Boundary Faces
    mesh.setMaxNumFaces( nBFa );
    mesh.numFaces() = nFa; // Here the REAL number of edges (all of them)
    mesh.setNumBFaces( nBFa );

    mesh.setMaxNumVolumes( nVo, true );

    mesh.setMarker( regionFlag ); // Add Marker to list of Markers

    typename RegionMesh3D<GeoShape, MC>::PointType * pp = 0;
    typename RegionMesh3D<GeoShape, MC>::EdgeType * pe = 0;
    typename RegionMesh3D<GeoShape, MC>::FaceType * pf = 0;
    typename RegionMesh3D<GeoShape, MC>::VolumeType * pv = 0;
    // addPoint()/Face()/Edge() returns a reference to the last stored point
    // I use that information to set all point info, by using a pointer.
    UInt count = 0;
    long int ibc;
    while ( next_good_line( mystream, line ).good() )
    {
        if ( line.find( "Vertices" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), mystream );
            for ( i = 0;i < nVe;i++ )
            {
                mystream >> x >> y >> z >> ibc;
                if ( ibc != 0 )
                {
                    ++count;
                    pp = &mesh.addPoint( true ); // Boundary point. Boundary switch set by the mesh method.
                    pp->setMarker( EntityFlag( ibc ) );
                }
                else
                {
                    pp = &mesh.addPoint( false );
                }
                pp->x() = x;
                pp->y() = y;
                pp->z() = z;
                pp->setMarker( EntityFlag( ibc ) );
            }
            oStr << "Vertices Read " << std::endl;
            done++;
            if ( count != nBVe )
                std::cerr << "NumB points inconsistent !" << std::endl;
        }

        if ( line.find( "Triangles" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), mystream );
            oStr << "Reading Bfaces " << std::endl;
            for ( i = 0;i < nBFa;i++ )
            {
                mystream >> p1 >> p2 >> p3 >> ibc;

                pf = &( mesh.addFace( true ) ); // Only boundary faces

                pf->setMarker( EntityFlag( ibc ) );
                pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.
            }
            oStr << "Boundary Faces Read " << std::endl;
            done++;
        }

        if ( line.find( "Quadrilaterals" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), mystream );
            oStr << "Reading Bfaces " << std::endl;
            for ( i = 0;i < nBFa;i++ )
            {
                mystream >> p1 >> p2 >> p3 >> p4 >> ibc;

                pf = &( mesh.addFace( true ) ); // Only boundary faces

                pf->setMarker( EntityFlag( ibc ) );
                pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.
                pf->setPoint( 4, mesh.point( p4 ) ); // set face conn.
            }
            oStr << "Boundary Faces Read " << std::endl;
            done++;
        }

        if ( line.find( "Edges" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), mystream );
            oStr << "Reading Bedges " << std::endl;
            for ( i = 0;i < nBEd;i++ )
            {
                mystream >> p1 >> p2 >> ibc;
                pe = &mesh.addEdge( true ); // Only boundary edges.
                pe->setMarker( EntityFlag( ibc ) );
                pe->setPoint( 1, mesh.point( p1 ) ); // set edge conn.
                pe->setPoint( 2, mesh.point( p2 ) ); // set edge conn.
            }
            oStr << "Boundary Edges Read " << std::endl;
            done++;
        }
        if ( line.find( "Tetrahedra" ) != std::string::npos )
        {
            count = 0;
            nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), mystream );
            oStr << "Reading Volumes " << std::endl;
            for ( i = 0; i < nVo; i++ )
            {
                mystream >> p1 >> p2 >> p3 >> p4 >> ibc;
                pv = &mesh.addVolume();
                pv->id() = i + 1;
                pv->setPoint( 1, mesh.point( p1 ) );
                pv->setPoint( 2, mesh.point( p2 ) );
                pv->setPoint( 3, mesh.point( p3 ) );
                pv->setPoint( 4, mesh.point( p4 ) );
                pv->setMarker( EntityFlag( ibc ) );
                count++;
            }
            oStr << count << " Volume elements Read" << std::endl;
            done++;
        }
        if ( line.find( "Hexahedra" ) != std::string::npos )
        {
            count = 0;
            nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), mystream );
            oStr << "Reading Volumes " << std::endl;
            for ( i = 0; i < nVo; i++ )
            {
                mystream >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> p8 >> ibc;
                pv = &mesh.addVolume();
                pv->id() = i + 1;
                pv->setPoint( 1, mesh.point( p1 ) );
                pv->setPoint( 2, mesh.point( p2 ) );
                pv->setPoint( 3, mesh.point( p3 ) );
                pv->setPoint( 4, mesh.point( p4 ) );
                pv->setPoint( 5, mesh.point( p5 ) );
                pv->setPoint( 6, mesh.point( p6 ) );
                pv->setPoint( 7, mesh.point( p7 ) );
                pv->setPoint( 8, mesh.point( p8 ) );
                pv->setMarker( EntityFlag( ibc ) );
                count++;
            }
            oStr << count << " Volume elements Read" << std::endl;
            done++;
        }

    }

    // Test mesh
    Switch sw;

    if ( !checkMesh3D( mesh, sw, true, verbose, oStr, std::cerr, oStr ) )
        abort();
    // if(!checkMesh3D(mesh, sw, true,true, oStr,oStr,oStr)) abort();//verbose version

    // This part is to build a P2 mesh from a P1 geometry

    if ( shape == TETRA && GeoShape::numPoints > 4 )
        p1top2( mesh );
    mystream.close();

    Real vols[ 3 ];
    getVolumeFromFaces( mesh, vols, oStr );
    oStr << "   VOLUME ENCLOSED BY THE MESH COMPUTED BY INTEGRATION ON" <<
    " BOUNDARY FACES" << std::endl;
    oStr << "INT(X)     INT(Y)      INT(Z) <- they should be equal and equal to" << std::endl <<
    "                                 the voulume enclosed by the mesh " << std::endl;
    oStr << vols[ 0 ] << " " << vols[ 1 ] << " " << vols[ 2 ] << std::endl;

    oStr << "   BOUNDARY FACES ARE DEFINING A CLOSED SURFACE IF " << testClosedDomain( mesh, oStr ) << std::endl <<
    " IS (ALMOST) ZERO" << std::endl;

    return done == 4 ;

}

//----------------------------------------------------------------------
//
// Problem: the functions that follows should be in a readMesh3D.cc
// nevertheless, if I do that, I can compile the library, but I
// obtain several errors during the linking with a main.cc.
// JFG, 07/09/2002
//======================================================================

// ****************      Mesh ++  Readers   **********************************
bool
readMppFileHead( std::ifstream & mystream,
                 UInt & numVertices,
                 UInt & numBVertices,
                 UInt & numBFaces,
                 UInt & numBEdges,
                 UInt & numVolumes );

// ****************      INRIA mesh  readers   **********************************

int nextIntINRIAMeshField( std::string const & line, std::istream & mystream );

//! Reads all basic info from INRIA MESH file
//! so as to be able to properly dimension all arrays
bool
readINRIAMeshFileHead( std::ifstream & mystream,
                       UInt & numVertices,
                       UInt & numBVertices,
                       UInt & numBFaces,
                       UInt & numBEdges,
                       UInt & numVolumes,
                       ReferenceShapes & shape );

}
#endif
