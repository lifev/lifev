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

#include <sstream>
#include <algorithm>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>

#include <debug.hpp>


#include "regionMesh3D.hpp"
#include "util_string.hpp"
#include "mesh_util.hpp"

namespace LifeV
{
/*----------------------------------------------------------------------*
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
readMppFileHead( std::ifstream & mystream,
                 UInt & numVertices,
                 UInt & numBVertices,
                 UInt & numBFaces,
                 UInt & numBEdges,
                 UInt & numVolumes );

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
        std::cerr << " Error in readMpp: File " << filename
                  << " not found or locked" << std::endl;
        abort();
    }
    std::cout << "Reading Mesh++ file" << std::endl;
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
        std::cerr << " Error in readMpp: File " << filename
                  << " not found or locked" << std::endl;
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
    oStr      << "#Faces    = "          << std::setw(10) << nFa
              << "  #Boundary Faces  = " << std::setw(10) << nBFa << std::endl;
    oStr      << "#Edges    = "          << std::setw(10) << nEd
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

/*!
  read an INRIA mesh
 */
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
        std::cerr << " Error in readINRIAMeshFile: File " << filename
                  << " not found or locked"
                  << std::endl;
        abort();
    }
    std::cout << "Reading INRIA mesh file" << std::endl;
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
        std::cerr << " Error in readINRIAMeshFile: File " << filename
                  << " not found or locked" << std::endl;
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
        std::cout << "Linear Hexa Mesh" << std::endl;
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
    oStr      << "#Faces    = "          << std::setw(10) << nFa
              << "  #Boundary Faces  = " << std::setw(10) << nBFa << std::endl;
    oStr      << "#Edges    = "          << std::setw(10) << nEd
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
//
// GMSH
//

template <typename GeoShape, typename MC>
bool
readGmshFile( RegionMesh3D<GeoShape, MC> & mesh,
             const std::string & filename,
             EntityFlag regionFlag )
{
    std::ifstream __is ( filename.c_str() );

    char __buf[256];
    __is >> __buf;
    Debug() << "buf: "<< __buf << "\n";
    uint __n;
    __is >> __n;
    Debug() << "number of nodes: " << __n;

    // Add Marker to list of Markers
    mesh.setMarker( regionFlag );


    std::vector<double> __x(3*__n);
    std::vector<bool> __isonboundary(__n);
    std::vector<uint> __whichboundary(__n);
    Debug() << "reading "<< __n << " nodes\n";
    std::map<int,int> itoii;
    for( uint __i = 0; __i < __n;++__i )
    {
		uint __ni;
		__is >> __ni
		     >> __x[3*__i]
		     >> __x[3*__i+1]
		     >> __x[3*__i+2];

        itoii[__ni-1] = __i;
    }
    __is >> __buf;
    Debug() << "buf: "<< __buf << "\n";
    __is >> __buf;
    Debug() << "buf: "<< __buf << "\n";
    uint __nele;
    __is >> __nele;

    typename RegionMesh3D<GeoShape, MC>::EdgeType * pe = 0;
    typename RegionMesh3D<GeoShape, MC>::FaceType * pf = 0;
    typename RegionMesh3D<GeoShape, MC>::VolumeType * pv = 0;




    Debug() << "number of elements: " << __nele << "\n";
    std::vector<std::vector<int> > __e(__nele);
    std::vector<int> __et(__nele);
    std::vector<int> __etype( __nele );
    std::vector<int> __gt(16);
    __gt.assign( 16, 0 );

    for( uint __i = 0; __i < __nele;++__i )
    {
		int __ne, __t, __tag, __np, __dummy;
		__is >> __ne
		     >> __t
		     >> __tag
		     >> __dummy
		     >> __np;


		++__gt[ __t ];
		__etype[__i] = __t;
		__et[__i] = __tag;
		__e[__i].resize( __np );
		int __p = 0;
		while ( __p != __np )
		{
		    __is >> __e[__i][__p];
            __e[__i][__p] = itoii[ __e[__i][__p]-1];
            __e[__i][__p] += 1;

		    ++__p;
		}
    }
    std::for_each( __gt.begin(), __gt.end(),  std::cout << boost::lambda::_1 << " " );
    std::cout << "\n";

    // Euler formulas
    uint n_volumes = __gt[4];
    uint n_faces_boundary = __gt[2];
    uint n_faces_total = 2*n_volumes+(n_faces_boundary/2);

    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges( __gt[1] );
    mesh.numEdges() = __gt[1]; // Here the REAL number of edges (all of them)
    mesh.setNumBEdges( __gt[1] );

    Debug() << "number of edges= " << __gt[1] << "\n";
    // Only Boundary Faces
    mesh.setMaxNumFaces( n_faces_total );
    mesh.numFaces() = n_faces_total; // Here the REAL number of edges (all of them)
    mesh.setNumBFaces( n_faces_boundary );
    Debug() << "number of faces= " << n_faces_total << "\n";
    mesh.setMaxNumVolumes( n_volumes, true );
    Debug() << "number of volumes= " << n_volumes << "\n";

    __isonboundary.assign( __n, false );
    __whichboundary.assign( __n, 0 );
    for( uint __i = 0; __i < __nele;++__i )
    {
        switch( __etype[__i] )
        {
            // triangular faces (linear)
            case 2:
            {
                __isonboundary[ __e[__i][0]-1 ] = true;
                __isonboundary[ __e[__i][1]-1 ] = true;
                __isonboundary[ __e[__i][2]-1 ] = true;

                __whichboundary[__e[__i][0]-1 ] = __et[__i];
                __whichboundary[__e[__i][1]-1 ] = __et[__i];
                __whichboundary[__e[__i][2]-1 ] = __et[__i];
            }
        }
    }
    // add the point to the mesh
    typename RegionMesh3D<GeoShape, MC>::PointType * pp = 0;

    mesh.setMaxNumPoints( __n, true );
    mesh.numVertices() = __n;
    mesh.numBVertices() = std::count( __isonboundary.begin(), __isonboundary.end(), true );
    mesh.setNumBPoints( mesh.numBVertices() );

    Debug() << "number of points : " << mesh.numPoints() << "\n";
    Debug() << "number of boundary points : " << mesh.numBPoints() << "\n";
    Debug() << "number of vertices : " << mesh.numVertices() << "\n";
    Debug() << "number of boundary vertices : " << mesh.numBVertices() << "\n";

    for( uint __i = 0; __i < __n;++__i )
    {
        pp = &mesh.addPoint( __isonboundary[ __i ] );
        pp->setMarker( __whichboundary[__i] );
        pp->x() = __x[3*__i];
        pp->y() = __x[3*__i+1];
        pp->z() = __x[3*__i+2];
    }

    // add the element to the mesh
    for( uint __i = 0; __i < __nele;++__i )
    {
        switch( __etype[__i] )
        {
            // segment(linear)
            case 1:
            {
                pe = &( mesh.addEdge( true ) );
                pe->setMarker( EntityFlag( __et[__i] ) );
                pe->setPoint( 1, mesh.point( __e[__i][0] ) );
                pe->setPoint( 2, mesh.point( __e[__i][1] ) );



            }
            break;
            // triangular faces (linear)
            case 2:
            {
                pf = &( mesh.addFace( true ) );
                pf->setMarker( EntityFlag( __et[__i] ) );
                pf->setPoint( 1, mesh.point( __e[__i][0] ) );
                pf->setPoint( 2, mesh.point( __e[__i][1] ) );
                pf->setPoint( 3, mesh.point( __e[__i][2] ) );

            }
            break;
            // quadrangular faces(linear)
            case 3:
            {
                pf = &( mesh.addFace( true ) );
                pf->setMarker( EntityFlag( __et[__i] ) );
                pf->setPoint( 1, mesh.point( __e[__i][0] ) );
                pf->setPoint( 2, mesh.point( __e[__i][1] ) );
                pf->setPoint( 3, mesh.point( __e[__i][2] ) );
                pf->setPoint( 4, mesh.point( __e[__i][3] ) );
            }
            break;
            // tetrahedrons(linear)
            case 4:
            {
                pv = &( mesh.addVolume() );
                pv->id() = __i + 1;
                pv->setMarker( EntityFlag( __et[__i] ) );
                pv->setPoint( 1, mesh.point( __e[__i][0] ) );
                pv->setPoint( 2, mesh.point( __e[__i][1] ) );
                pv->setPoint( 3, mesh.point( __e[__i][2] ) );
                pv->setPoint( 4, mesh.point( __e[__i][3] ) );
            }
            break;
            // hexahedrons(linear)
            case 5:
            {
                pv = &( mesh.addVolume() );

                pv->id() = __i + 1;
                pv->setMarker( EntityFlag( __et[__i] ) );
                pv->setPoint( 1, mesh.point( __e[__i][0] ) );
                pv->setPoint( 2, mesh.point( __e[__i][1] ) );
                pv->setPoint( 3, mesh.point( __e[__i][2] ) );
                pv->setPoint( 4, mesh.point( __e[__i][3] ) );
                pv->setPoint( 5, mesh.point( __e[__i][4] ) );
                pv->setPoint( 6, mesh.point( __e[__i][5] ) );
                pv->setPoint( 7, mesh.point( __e[__i][6] ) );
                pv->setPoint( 8, mesh.point( __e[__i][7] ) );
            }
            break;

        }


    }


    Switch sw;
    if( checkMesh3D(mesh, sw, true,true,std::cout,std::cout,std::cout) == false )
    {
        std::ostringstream __ex;
        __ex << "invalid mesh from GSMH";
        throw std::logic_error( __ex.str() );
    }

    return true;
}



template<typename GeoShape, typename MC>
bool
readNetgenMesh(RegionMesh3D<GeoShape,MC> & mesh,
               const std::string  & filename,
               EntityFlag regionFlag);





/*
  MM added support for using meshes
  generated by netgen

  I referred a bit to: 
  LifeV/.../mesh/readMesh3D.h/.cc
  netgen/libsrc/meshing/meshclass.cpp

#include "markers.h"

#include "regionMesh3D.h"
#include "util_string.h"
#include "mesh_util.h"
*/

#include "bareItems.hpp"

template<typename GeoShape, typename MC>
bool
readNetgenMesh(RegionMesh3D<GeoShape,MC> & mesh, 
               const std::string  & filename,
               EntityFlag regionFlag)
{
  std::string line;
  UInt nVe(0), nBVe(0),nFa(0),nBFa(0),nPo(0),nBPo(0);
  UInt nVo(0),nEd(0),nBEd(0);
  UInt i;
  UInt fake;
  UInt surfnr,bcnr,domin,domout;
  UInt surf1,surf2,matnr,np,p1,p2,p3,p4,t;
  Real x,y,z;
  std::vector<bool> bpoints;
//  std::vector<bool> bedges;
  BareItemsHandler<BareEdge> bihBedges;
  std::vector<EntityFlag> bcnsurf,bcnpoints;
  UInt flag;
  typename MC::PointMarker PMarker;
  typename MC::EdgeMarker EMarker;


  std::ifstream fstream1(filename.c_str()); 
  if (fstream1.fail()) {
    std::cerr<<"Error in readNetgenMesh: File not found or locked"<<std::endl;
    abort();
  }
  std::cout<< "Reading netgen mesh file: "<<filename<<std::endl;
  fstream1>>line;
  if(line!="mesh3d"){
    std::cerr<<"Error in readNetgenMesh: mesh file is not in mesh3d format (netgen)"\
        <<std::endl;
    abort();
  }

  /* I assume as I tested that faces stored are only boundary faces 
     and edges stored are only a part of boundary ones with
     inside a face with a common marker, instead points can be
     inside the domain too, but this format doesn't say me
     which, so I'll find them from myself
   */
  flag=1|2|4|8;
  while(!fstream1.eof())
  {
    fstream1>>line;

    if(line=="points" && flag&1){
      fstream1>>nVe;
      bcnpoints.reserve(nVe+1);
      bpoints.reserve(nVe+1);
      for(i=0;i<nVe+1;i++){
        bpoints[i]=false;
        bcnpoints[i]=NULLFLAG;
      }
      flag&=~1;
      break;
    }
  }
  fstream1.close();
  
/* as Forma told, I have found sometimes problems with seekg, becouse
   it seems that sometimes it put the stream in a non consistent state,
   so I open new files instead of seeking, I have had not compile errors
   only sometimes problems on big mesh files with seekg
*/
//  fstream.seekg(0,std::ios_base::beg);
  std::ifstream fstream2(filename.c_str()); 
  while(!fstream2.eof())
  {
    fstream2>>line;

    if(line=="surfaceelementsgi" && flag&8){
      fstream2>>nBFa;
      bcnsurf.reserve(nBFa+1);
      bcnsurf[0]=NULLFLAG;
      for(i=0;i<nBFa;i++){
        fstream2>>surfnr>>bcnr>>domin>>domout>>np>>p1>>p2>>p3\
               >>t>>t>>t;
        bcnsurf[i+1]=bcnr;
        if(np!=3){
          std::cerr<<"Error in readNetgenMesh: " \
                "only triangular surfaces supported"<<std::endl;
          abort();
        }
        //assume p1!=p2!=p3!=p1
        nBVe+=bpoints[p1]?0:1;
        nBVe+=bpoints[p2]?0:1;
        nBVe+=bpoints[p3]?0:1;
        bpoints[p1]=bpoints[p2]=bpoints[p3]=true;

        /* here I set the boundary points marker
           note: this works only with my patch
           of strongerFlag */
        bcnpoints[p1]=PMarker.setStrongerMarker(bcnpoints[p1],bcnr);
        bcnpoints[p2]=PMarker.setStrongerMarker(bcnpoints[p2],bcnr);
        bcnpoints[p3]=PMarker.setStrongerMarker(bcnpoints[p3],bcnr);

/* now I have the surface and points, so I can calculate
   the num of edges on the boundary, useful in case this is
   a quadratic tetra, I calculate the bcnr too using
   BareItemHandler<BareEdge> to store results
*/
   //I've found this silly but easy way
        BareEdge bed=setBareEdge(p1,p2);
        bihBedges.addIfNotThere(bed,(ID)NULLFLAG);
        bihBedges[bed]=(ID)EMarker.setStrongerMarker(
                              bihBedges[bed],
                              bcnr);
        bed=setBareEdge(p2,p3);
        bihBedges.addIfNotThere(bed,(ID)NULLFLAG);
        bihBedges[bed]=(ID)EMarker.setStrongerMarker(
                              bihBedges[bed],
                              bcnr);
        bed=setBareEdge(p3,p1);
        bihBedges.addIfNotThere(bed,(ID)NULLFLAG);
        bihBedges[bed]=(ID)EMarker.setStrongerMarker(
                              bihBedges[bed],
                              bcnr);
      }
      flag&=~8;
      break;
    }
  }
  fstream2.close();
  nBEd=bihBedges.howMany();
  


/* I assume as I tested that edges stored are only a part of
   boundary ones, so really they are meaningless to me */
  std::ifstream fstream3(filename.c_str()); 
  while(!fstream3.eof())
  {
    fstream3>>line;
    if(line=="edgesegmentsgi2" && flag&2){
      fstream3>>fake;		//this is not really the nBEd
      fake*=2;		//there are twice the num of lines told, test if prob
      for(i=0;i<fake;i++){
        fstream3>>surf1>>surf2>>p1>>p2>>t>>t>>t>>t>>t>>t>>t>>t;
      }
      flag&=~2;
    }
    else if(line=="volumeelements" && flag&4){
      fstream3>>nVo;
      for(i=0;i<nVo;i++){
        fstream3>>matnr>>np>>p1>>p2>>p3>>p4;
        if(np!=4){
          std::cerr<<"Error in readNetgenMesh: " \
                "only tetrahedra elements supported"<<std::endl;
          abort();
        }
      }
      flag&=~4;
    }
    if(!flag)break;
  }
  fstream3.close();
  if(flag!=0){
    std::cerr<<"Error in readNetgenMesh: " \
          "the mesh file does not have all the required sections" \
        <<std::endl;
    abort();
  }

  // Euler formulas
  nFa=2*nVo+(nBFa/2);
  nEd=nVo+nVe+(3*nBFa-2*nBVe)/4;

  // Be a little verbose
  if (GeoShape::numPoints > 4 ){
    std::cout << "Quadratic Tetra  Mesh (from Linear geometry)" <<std::endl;
    nPo=nVe+nEd;
    nBPo=nBVe+nBEd;	// I calculated the real nBEd before
  } else {
    std::cout << "Linear Tetra Mesh" <<std::endl;
    nPo=nVe;
    nBPo=nBVe;		
  }

  //points can be only vertices or on edges too

  std::cout<< "#Vertices= "<<nVe;
  std::cout<< " #BVertices= "<<nBVe<<std::endl;
  std::cout<< "#Faces= "<<nFa;
  std::cout<< " #Boundary Faces= "<<nBFa<<std::endl;
  std::cout<< "#Edges= "<<nEd;
  std::cout<< " #Boundary Edges= "<<nBEd<<std::endl;
  std::cout<< "#Points= "<<nPo;
  std::cout<< " #Boundary Points= "<<nBPo<<std::endl;
  std::cout<< "#Volumes= "<<nVo<<std::endl;

  // Set all basic data structure

  // I store all Points
  mesh.setMaxNumPoints(nPo,true);
  mesh.setNumBPoints(nBPo);
  mesh.numVertices()=nVe;
  mesh.numBVertices()=nBVe;
  // Only Boundary Edges (in a next version I will allow for different choices)
  mesh.setMaxNumEdges(nBEd);
  mesh.numEdges()=nEd; // Here the REAL number of edges (all of them)
  mesh.setNumBEdges(nBEd);	/////////????????????
  // Only Boundary Faces
  mesh.setMaxNumFaces(nBFa);
  mesh.numFaces()=nFa; // Here the REAL number of edges (all of them)
  mesh.setNumBFaces(nBFa);

  mesh.setMaxNumVolumes(nVo,true);

  mesh.setMarker(regionFlag); // Mark the region ????????what if more then one<<<<<<<

  typename RegionMesh3D<GeoShape,MC>::PointType * pp=0;
  typename RegionMesh3D<GeoShape,MC>::EdgeType * pe=0;
  typename RegionMesh3D<GeoShape,MC>::FaceType * pf=0;
  typename RegionMesh3D<GeoShape,MC>::VolumeType * pv=0;
  // addPoint()/Face()/Edge() returns a reference to the last stored point
  // I use that information to set all point info, by using a pointer.


  std::ifstream  fstream4(filename.c_str()); 
  flag=1|2|4|8;
  while(!fstream4.eof())
  {
    fstream4>>line;

    if(line=="points" && flag&1){
      fstream4>>nVe;
      for(i=0;i<nVe;i++){
        fstream4>>x>>y>>z;
        pp=&mesh.addPoint(bpoints[i+1]); //true if boundary point
        pp->setMarker(bcnpoints[i+1]); 
        pp->x()=x;
        pp->y()=y;
        pp->z()=z;       
      }
      flag&=~1;
      break;
    }
  }
  fstream4.close();
  std::ifstream  fstream5(filename.c_str()); 

  while(flag && !fstream5.eof())  //fstream after using .seekg sometimes never reaches eof, 
  				//entering infinite loop, why????? damn STL I reopen the file and put the flag to fix
  {
    fstream5>>line;

    if(line=="edgesegmentsgi2" && flag&2){
      fstream5>>fake;	//this is not really the nBEd
      fake*=2;		//there are twice the num of lines told, test if prob
      for(i=0;i<fake;i++){
        fstream5>>surf1>>surf2>>p1>>p2>>t>>t>>t>>t>>t>>t>>t>>t;
      }
      /* here I set the real boundary edges that I stored 
         in bihBedges
       */
      BareItemsHandler<BareEdge>::const_iterator bedge=bihBedges.begin();
      for(i=0;i<nBEd;i++){
        pe = &mesh.addEdge( true ); // Only boundary edges.
        pe->setMarker( EntityFlag(bedge->second)); 
        p1=bedge->first.first;
        p2=bedge->first.second;
        pe->setPoint( 1, mesh.point( p1 ) ); // set edge conn.
        pe->setPoint( 2, mesh.point( p2 ) ); // set edge conn.
        bedge++;
      }
      flag&=~2;
    }
    else if(line=="volumeelements" && flag&4){
      fstream5>>nVo;
      for(i=0;i<nVo;i++){
        fstream5>>matnr>>np>>p1>>p2>>p3>>p4;
        pv=&mesh.addVolume();
        pv->id()=i+1;
        pv->setPoint(1, mesh.point(p1) );
        pv->setPoint(2, mesh.point(p2) );
        pv->setPoint(3, mesh.point(p3) );
        pv->setPoint(4, mesh.point(p4) );
      }
      flag&=~4;
    }
    else if(line=="surfaceelementsgi" && flag&8){
      fstream5>>nFa;
      for(i=0;i<nFa;i++){
        fstream5>>surfnr>>bcnr>>domin>>domout>>np>>p1>>p2>>p3\
               >>t>>t>>t;
        pf=&mesh.addFace(true); // Only boundary faces

        pf->setMarker(EntityFlag(bcnr));
        pf->setPoint(1,mesh.point(p1)); // set face conn.
        pf->setPoint(2,mesh.point(p2)); // set face conn.
        pf->setPoint(3,mesh.point(p3)); // set face conn.
      }
      flag&=~8;
    }
  }

  // This part is to build a P2 mesh from a P1 geometry
  
  if (GeoShape::numPoints > 4 )p1top2(mesh);
  fstream5.close();
  
  // Test mesh
  Switch sw;
  
  ///// CORRECTION JFG
  //if (mesh.check(1, true,true))done=0;
  
  if(!checkMesh3D(mesh, sw, true,true,std::cout,std::cout,std::cout)) abort(); // CORRECTION JFG

  Real vols[3];
  getVolumeFromFaces(mesh, vols,std::cout);
  std::cout<< "   VOLUME ENCLOSED BY THE MESH COMPUTED BY INTEGRATION ON"<<
    " BOUNDARY FACES"<<std::endl;
  std::cout << "INT(X)     INT(Y)      INT(Z) <- they should be equal and equal to"<<std::endl<<
    "                                 the voulume enclosed by the mesh "<<std::endl;
  std::cout<<vols[0]<<" "<<vols[1]<<" "<<vols[2]<<std::endl;
  
  std::cout<< "   BOUNDARY FACES ARE DEFINING A CLOSED SURFACE IF "<<testClosedDomain(mesh,std::cout)<< std::endl<<
    " IS (ALMOST) ZERO"<<std::endl;
  
  return true;
}


#include <iostream>
#include <fstream>

/* 
MM I've looked to this source, easy func
/usr/local/src/ng431/libsrc/interface/importsolution.cpp
*/
template <typename VectorType>
void saveNetgenSolution(std::string filename,const VectorType& U,std::string fctname="u")
{
  std::ofstream of(filename.c_str());
  of<<"solution "<<fctname<<" -size="<<U.size()
    <<" -components=1 -type=nodal"<<std::endl;
  for(UInt i=0;i<U.size();i++)
    of<<U(i)<<std::endl;
  of.close();
}

}
#endif
