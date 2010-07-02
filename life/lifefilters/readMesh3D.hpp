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
#ifndef _READMESH3D_HH_
#define _READMESH3D_HH_

#include <sstream>
#include <algorithm>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>

#include <life/lifecore/debug.hpp>


#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifecore/util_string.hpp>
#include <life/lifefilters/mesh_util.hpp>
#include <life/lifefilters/selectMarker.hpp>

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

// template<typename T>
// T Max(T & i1, T & i2,T & i3)
// {
//     return std::max(std::max(i1,i2),i3);
// }

// template<typename T>
// T Max(T & i1, T & i2, T & i3, T & i4)
// {
//     return std::max(Max(i1,i2,i3),i4);
// }

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
             bool verbose = false )
{
    unsigned done = 0;
    std::string line;
    Real x, y, z;
    int ity, ity_id;
    UInt p1, p2, p3, p4;
    UInt nVe( 0 ), nBVe( 0 ), nFa( 0 ), nBFa( 0 ), nPo( 0 ), nBPo( 0 );
    UInt nVo( 0 ), nBEd( 0 );
    UInt nEd( 0 );
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
    mesh.setMaxNumGlobalPoints( nPo );
    mesh.setNumBPoints  ( nBPo );
    mesh.setNumVertices ( nVe );
    mesh.setNumGlobalVertices(nVe);
    mesh.setNumBVertices( nBVe );
    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges ( nEd );
    mesh.setMaxNumGlobalEdges ( nEd );
    mesh.setNumEdges    ( nEd ); // Here the REAL number of edges (all of them)
    mesh.setNumBEdges   ( nBEd );
    // Only Boundary Faces
    mesh.setMaxNumFaces ( nBFa );
    mesh.setMaxNumGlobalFaces ( nBFa );
    mesh.setNumFaces    ( nFa ); // Here the REAL number of edges (all of them)
    mesh.setNumBFaces   ( nBFa );

    mesh.setMaxNumVolumes( nVo, true );
    mesh.setMaxNumGlobalVolumes( nVo);

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

                pp->setId     ( i + 1 );
                pp->setLocalId( i + 1 );

                mesh.localToGlobalNode().insert(std::make_pair(i+1, i+1));
                mesh.globalToLocalNode().insert(std::make_pair(i+1, i+1));
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
                pv->setId     ( i + 1 );
                pv->setLocalId( i + 1);
//                pv->id() = i + 1;
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
                       UInt & numStoredFaces,
                       ReferenceShapes & shape,
                       InternalEntitySelector
                       iSelect=InternalEntitySelector());

/*!
  read an INRIA mesh
*/
struct FiveNumbers
{
public:
    UInt i1,i2,i3,i4;
    long int ibc;
};
template <typename GeoShape, typename MC>
bool
readINRIAMeshFile( RegionMesh3D<GeoShape, MC>&      mesh,
                   std::string const&               filename,
                   EntityFlag                       regionFlag,
                   bool verbose                   = false,
                   InternalEntitySelector iSelect = InternalEntitySelector())
{
    unsigned done = 0;
    std::string line;
    Real x, y, z;
    UInt p1, p2, p3, p4, p5, p6, p7, p8;
    UInt nVe( 0 ), nFa( 0 ), nBFa( 0 ), nPo( 0 ), nBPo( 0 );
    UInt nBVe(0);
    UInt numStoredFaces(0);
    UInt nVo( 0 ), nBEd( 0 );
    UInt nEd;
    UInt i;
    ReferenceShapes shape(NONE);
    std::vector<FiveNumbers> faceHelp;
    typename std::vector<FiveNumbers>::iterator faceHelpIterator;
    std::stringstream discardedLog;
    std::ostream& oStr = verbose ? std::cout : discardedLog;

    // open stream to read header

    std::ifstream hstream( filename.c_str() );
    if (verbose)
        {
            std::cout<<"Reading form file "<<filename<< std::endl;
        }

    if ( hstream.fail() )
    {
        std::cerr << " Error in readINRIAMeshFile: File " << filename
                  << " not found or locked"
                  << std::endl;
        abort();
    }
    if (verbose) std::cout << "Reading INRIA mesh file" << filename << std::endl;
    if ( ! readINRIAMeshFileHead( hstream, nVe, nBVe, nBFa, nBEd, nVo, numStoredFaces,shape,
                                  iSelect) )
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
    int num1  = nVe + nVo;
    int num2  = nBVe;
    int num3  = nBFa;

    nEd = (3*num3 - 2*num2)/4 + num1;

//    nEd = (int) nVo + nVe + ( 3 * nBFa + dummy ) / 4;

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
                nBEd=(int( nBVe + nBFa - int(2) )>0?( nBVe + nBFa - 2 ):0);
                nBPo = (int(nBVe + ( nBVe + nBFa - int(2) ))>0?nBVe + ( nBVe + nBFa - 2 ):0);
            }
            else
            {
                if (verbose)
                    std::cout << "Linear Tetra Mesh" << std::endl;

                nPo = nVe;
                nBPo = nBVe;
                nBEd=(int( nBVe + nBFa - int(2) )>0?( nBVe + nBFa - 2 ):0);
            }
            break;
        default:
            ERROR_MSG( "Current version of INRIA Mesh file reader only accepts TETRA and HEXA" );
    }

    oStr << "#Vertices = "          << std::setw(10) << nVe
              << "  #BVertices       = " << std::setw(10) << nBVe << std::endl;
    oStr << "#Faces    = "          << std::setw(10) << nFa
              << "  #Boundary Faces  = " << std::setw(10) << nBFa << std::endl
              << "#Stored Faces = " << std::setw(10) << numStoredFaces<< std::endl;
    oStr << "#Edges    = "          << std::setw(10) << nEd
              << "  #Boundary Edges  = " << std::setw(10) << nBEd << std::endl;
    oStr << "#Points   = "          << std::setw(10) << nPo
              << "  #Boundary Points = " << std::setw(10) << nBPo << std::endl;
    oStr << "#Volumes  = "          << std::setw(10) << nVo  << std::endl;

    // Set all basic data structure

    // I store all Points
    mesh.setMaxNumPoints   ( nPo, true );
    mesh.setMaxNumGlobalPoints( nPo );
    mesh.setNumBPoints     ( nBPo );
    mesh.setNumVertices    ( nVe );
    mesh.setNumGlobalVertices(nVe);
    mesh.setNumBVertices   ( nBVe );
    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges    ( nBEd );
    mesh.setMaxNumGlobalEdges ( nEd );
    mesh.setNumEdges       ( nEd ); // Here the REAL number of edges (all of them)
    mesh.setNumBEdges      ( nBEd );
    // Only Boundary Faces
    mesh.setMaxNumFaces    ( numStoredFaces );
    mesh.setMaxNumGlobalFaces ( nBFa );
    mesh.setNumFaces       ( nFa ); // Here the REAL number of faces (all of them)
    mesh.setNumBFaces      ( nBFa );

    mesh.setMaxNumVolumes  ( nVo, true );
    mesh.setMaxNumGlobalVolumes( nVo);

    mesh.setMarker         ( regionFlag ); // Add Marker to list of Markers

    typedef typename RegionMesh3D<GeoShape, MC>::PointType  PointType;
    typedef typename RegionMesh3D<GeoShape, MC>::VolumeType VolumeType;


    typename RegionMesh3D<GeoShape, MC>::PointType * pp = 0;
    typename RegionMesh3D<GeoShape, MC>::EdgeType * pe = 0;
    typename RegionMesh3D<GeoShape, MC>::FaceType * pf = 0;
    typename RegionMesh3D<GeoShape, MC>::VolumeType * pv = 0;
    // addPoint()/Face()/Edge() returns a reference to the last stored point
    // I use that information to set all point info, by using a pointer.
    UInt count = 0;
    long int ibc;

    // To account for internal faces
    if (numStoredFaces > nBFa){
        faceHelp.resize(numStoredFaces - nBFa);
        faceHelpIterator=faceHelp.begin();
        oStr<<"WARNING: The mesh file (apparently) contains "<<numStoredFaces - nBFa<<" internal faces"<<std::endl;

    }

    while ( next_good_line( mystream, line ).good() )
    {
        if ( line.find( "Vertices" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), mystream );
            for ( i = 0; i < nVe; i++ )
            {
                mystream >> x >> y >> z >> ibc;

//                if (ibc == 1 ) ibc = 100;

                if ( !iSelect(EntityFlag(ibc)))
                {
                    ++count;
                    pp = &mesh.addPoint( true ); // Boundary point. Boundary switch set by the mesh method.
                    pp->setMarker( EntityFlag( ibc ) );
                }
                else
                {
                    pp = &mesh.addPoint( false );
                }
                pp->setId     ( i + 1 );
                pp->setLocalId( i + 1 );
                pp->x() = x;
                pp->y() = y;
                pp->z() = z;
                pp->setMarker( EntityFlag( ibc ) );

                mesh.localToGlobalNode().insert(std::make_pair(i+1, i+1));
                mesh.globalToLocalNode().insert(std::make_pair(i+1, i+1));
            }
            oStr << "Vertices Read " << std::endl;
            oStr << "size of the node storage is " << count*sizeof(PointType)/1024./1024. << std::endl;
            done++;
            if ( count != nBVe )
                std::cerr << "NumB points inconsistent !" << std::endl;
        }

        if ( line.find( "Triangles" ) != std::string::npos ){
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), mystream );
            oStr << "Reading Bfaces " << std::endl;
            for ( i = 0;i < numStoredFaces;i++ )
            {
                mystream >> p1 >> p2 >> p3 >> ibc;

                if (numStoredFaces > nBFa){
                    if (mesh.point( p1 ).boundary()&&mesh.point( p2 ).boundary()&&
                        mesh.point( p3 ).boundary()){
                        pf = &( mesh.addFace( true ) ); // Boundary faces
                        pf->setMarker( EntityFlag( ibc ) );
                        pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                        pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                        pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.

                    } else {
                        faceHelpIterator->i1=p1;
                        faceHelpIterator->i2=p2;
                        faceHelpIterator->i3=p3;
                        faceHelpIterator->ibc=ibc;
                        ++faceHelpIterator;
                    }
                } else {

                    pf = &( mesh.addFace( true ) ); // Only boundary faces
                    pf->setMarker( EntityFlag( ibc ) );
                    pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                    pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                    pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.
                }
            }
            for (faceHelpIterator=faceHelp.begin();faceHelpIterator!=faceHelp.end();
                 ++faceHelpIterator){
                p1=faceHelpIterator->i1;
                p2=faceHelpIterator->i2;
                p3=faceHelpIterator->i3;
                ibc=faceHelpIterator->ibc;
                pf = &( mesh.addFace( false ) ); // INTERNAL FACE
                pf->setMarker( EntityFlag( ibc ) );
                pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.
            }

            oStr << "Boundary Faces Read " << std::endl;
            done++;
        }

        if ( line.find( "Quadrilaterals" ) != std::string::npos ){
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), mystream );
            oStr << "Reading Bfaces " << std::endl;
            for ( i = 0;i < nBFa;i++ )
            {
                mystream >> p1 >> p2 >> p3 >> p4 >> ibc;

                if (numStoredFaces > nBFa){
                    if (mesh.point( p1 ).boundary()&&mesh.point( p2 ).boundary()&&
                        mesh.point( p3 ).boundary()){
                        pf = &( mesh.addFace( true ) ); // Boundary faces
                        pf->setMarker( EntityFlag( ibc ) );
                        pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                        pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                        pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.
                        pf->setPoint( 4, mesh.point( p4 ) ); // set face conn.

                    } else {
                        faceHelpIterator->i1=p1;
                        faceHelpIterator->i2=p2;
                        faceHelpIterator->i3=p3;
                        faceHelpIterator->i4=p4;
                        faceHelpIterator->ibc=ibc;
                        ++faceHelpIterator;
                    }
                } else {
                    pf = &( mesh.addFace( true ) ); // Only boundary faces
                    pf->setMarker( EntityFlag( ibc ) );
                    pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                    pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                    pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.
                    pf->setPoint( 4, mesh.point( p4 ) ); // set face conn.
                }
            }
            oStr << "Boundary Faces Read " << std::endl;
            for (faceHelpIterator=faceHelp.begin();faceHelpIterator!=faceHelp.end();
                 ++faceHelpIterator){
                p1=faceHelpIterator->i1;
                p2=faceHelpIterator->i2;
                p3=faceHelpIterator->i3;
                p4=faceHelpIterator->i4;
                ibc=faceHelpIterator->ibc;
                pf = &( mesh.addFace( false ) ); // INTERNAL FACE
                pf->setMarker( EntityFlag( ibc ) );
                pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.
                pf->setPoint( 4, mesh.point( p4 ) ); // set face conn.
            }
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
                pv->setId     ( i + 1 );
                pv->setLocalId( i + 1);
                pv->setPoint( 1, mesh.point( p1 ) );
                pv->setPoint( 2, mesh.point( p2 ) );
                pv->setPoint( 3, mesh.point( p3 ) );
                pv->setPoint( 4, mesh.point( p4 ) );
                pv->setMarker( EntityFlag( ibc ) );
//                mesh.localToGlobalElem().insert(std::make_pair(i+1, i+1));
//                mesh.globalToLocalElem().insert(std::make_pair(i+1, i+1));
                count++;
            }
            oStr << "size of the volume storage is "
                 << sizeof(VolumeType)*count/1024./1024. << " Mo." << std::endl;
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
                pv->setId     ( i + 1 );
                pv->setLocalId( i + 1);
//                pv->id() = i + 1;
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

/**
   read a gmsh mesh (3D) file and store it in a RegionMesh3D
   @param mesh mesh data structure to fill in
   @param filename name of the gmsh mesh file  to read
   @param regionFlag identifier for the region
   @param verbose whether the function shall be verbose
   @return true if everything went fine, false otherwise
*/
template <typename GeoShape, typename MC>
bool
readGmshFile( RegionMesh3D<GeoShape, MC> & mesh,
              const std::string & filename,
              EntityFlag regionFlag,
              bool verbose=false )
{
    std::ifstream __is ( filename.c_str() );
    Debug() << "gmsh reading: "<< filename << "\n";

    //    char __buf[256];
    std::string __buf;
    for (int ii = 0; ii < 6; ++ii)
    {
        __is >> __buf;
        std::cout << "buf: "<< __buf << "\n";
    }

    UInt __n;
    __is >> __n;
    Debug() << "number of nodes: " << __n;

    // Add Marker to list of Markers
    mesh.setMarker( regionFlag );


    std::vector<double> __x(3*__n);
    std::vector<bool> __isonboundary(__n);
    std::vector<UInt> __whichboundary(__n);
    Debug() << "reading "<< __n << " nodes\n";
    std::map<int,int> itoii;

    for( UInt __i = 0; __i < __n;++__i )
    {
        UInt __ni;
        __is >> __ni
             >> __x[3*__i]
             >> __x[3*__i+1]
             >> __x[3*__i+2];

        itoii[__ni - 1] = __i;
    }
    __is >> __buf;
    Debug() << "buf: "<< __buf << "\n";
    __is >> __buf;
    Debug() << "buf: "<< __buf << "\n";
    UInt __nele;
    __is >> __nele;

    typename RegionMesh3D<GeoShape, MC>::EdgeType   * pe = 0;
    typename RegionMesh3D<GeoShape, MC>::FaceType   * pf = 0;
    typename RegionMesh3D<GeoShape, MC>::VolumeType * pv = 0;




    Debug() << "number of elements: " << __nele << "\n";
    std::vector<std::vector<int> > __e(__nele);
    std::vector<int> __et(__nele);
    std::vector<int> __etype( __nele );
    std::vector<int> __gt(32);
    __gt.assign( 32, 0 );

    for( UInt __i = 0; __i < __nele;++__i )
    {
        int __ne, __t, __np;

        //Debug() << __i + 1 << " ";

        __is >> __buf;
        __is >> __ne;

        //Debug() << __ne << " ";

        switch (__ne)
        {
        case(2):
            __np = 3;
            break;
        case(4):
            __np = 4;
            break;
        case(15):
            __np = 1;
            break;
        default:
            __np = 0;
            Debug() << "Element type unknown " << __ne << "\n";
            ASSERT(true, "Elements type unsupported.\n")
        }

        __is >> __t;

        //Debug() << __t << " ";

        bool ibcSet = false;
        int  flag   = 0;
        int __tag(0);

        for (int iflag = 0; iflag < __t; ++iflag)
        {
            __is >> flag;

            if (!ibcSet)
            {
                __tag = flag;
                ibcSet = true;
            }
        }

        ++__gt[ __ne];

        __etype[__i] = __ne;

        __et[__i] = __tag;
        __e[__i].resize( __np );

        int __p = 0;
        while ( __p != __np )
        {
            int node;
            __is >> node;
            __e[__i][__p] = node;
            __e[__i][__p] = itoii[ __e[__i][__p] - 1];
            __e[__i][__p] += 1;

            ++__p;
        }
    }


    // Euler formulas
    UInt n_volumes = __gt[4];
    UInt n_faces_boundary = __gt[2];
    UInt n_faces_total = 2*n_volumes+(n_faces_boundary/2);

    mesh.setMaxNumGlobalPoints( __n );
    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges( __gt[1] );
    mesh.setNumEdges   ( __gt[1] ); // Here the REAL number of edges (all of them)
    mesh.setNumBEdges  ( __gt[1] );
    mesh.setMaxNumGlobalEdges( __gt[1] );
    Debug() << "number of edges= " << __gt[1] << "\n";
    // Only Boundary Faces
    mesh.setMaxNumFaces( n_faces_total );
    mesh.setNumFaces   ( n_faces_total ); // Here the REAL number of edges (all of them)
    //mesh.setMaxNumFaces( n_faces_boundary );
    //mesh.setNumFaces   ( n_faces_boundary ); // Here the REAL number of edges (all of them)
    mesh.setNumBFaces  ( n_faces_boundary );
    mesh.setMaxNumGlobalFaces( n_faces_total );
    Debug() << "number of faces= " << n_faces_boundary << "\n";
    mesh.setMaxNumVolumes( n_volumes, true );
    mesh.setMaxNumGlobalVolumes( n_volumes );
    Debug() << "number of volumes= " << n_volumes << "\n";

    __isonboundary.assign( __n, false );
    __whichboundary.assign( __n, 0 );
    for( UInt __i = 0; __i < __nele;++__i )
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
    mesh.setNumVertices ( __n );
    mesh.setNumBVertices( std::count( __isonboundary.begin(), __isonboundary.end(), true ) );
    mesh.setNumBPoints  ( mesh.numBVertices() );

    Debug() << "number of points : " << mesh.numPoints() << "\n";
    Debug() << "number of boundary points : " << mesh.numBPoints() << "\n";
    Debug() << "number of vertices : " << mesh.numVertices() << "\n";
    Debug() << "number of boundary vertices : " << mesh.numBVertices() << "\n";

    for( UInt __i = 0; __i < __n;++__i )
    {
        pp = &mesh.addPoint( __isonboundary[ __i ] );
        pp->setMarker( __whichboundary[__i] );
        pp->setId     ( __i + 1 );
        pp->setLocalId( __i + 1 );
        pp->x() = __x[3*__i];
        pp->y() = __x[3*__i+1];
        pp->z() = __x[3*__i+2];
        mesh.localToGlobalNode().insert(std::make_pair(__i+1, __i+1));
        mesh.globalToLocalNode().insert(std::make_pair(__i+1, __i+1));
    }

    int nVo = 1;
    // add the element to the mesh
    for( UInt __i = 0; __i < __nele;++__i )
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
                pv->setId     ( nVo );
                pv->setLocalId( nVo++);
//                pv->id() = __i + 1;
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

                pv->setId     ( __i + 1 );
                pv->setLocalId( __i + 1 );
//                pv->id() = __i + 1;
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
    std::stringstream discardedLog;
    std::ostream& oStr = verbose ? std::cout : discardedLog;
    if( checkMesh3D(mesh, sw, true,verbose,oStr,std::cerr,oStr) == false )
    {
        std::ostringstream __ex;
        __ex << "invalid mesh from GSMH";
        throw std::logic_error( __ex.str() );
    }

    return true;
}



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

#include <life/lifemesh/bareItems.hpp>

template<typename GeoShape, typename MC>
bool
readNetgenMesh(RegionMesh3D<GeoShape,MC> & mesh,
               const std::string  & filename,
               EntityFlag regionFlag,
               bool verbose=false )
{
	// I will extract lines from iostream
    std::string line;
    // number of Geo Elements
    UInt nVe(0), nBVe(0), nPo(0), nBPo(0);
    UInt nEd(0), nBEd(0), nFa(0), nBFa(0);
    UInt nVo(0);
    // During the first access to file, build list of structures
    Vector pointcoor;
    std::vector<UInt> edgepointID, facepointID, volumepointID;
    // for poit i-th: is it a boindary point?
    std::vector<bool> bpoints;
    // build a list of boundary edges, since netgen is not writing all of them
    BareItemsHandler<BareEdge> bihBedges;
    // flags for boundary entities
    std::vector<EntityFlag> bcnsurf, bcnpoints;
    // bitstream to check which file section has already been visited
    UInt flag;

    typename MC::PointMarker PMarker;
    typename MC::EdgeMarker EMarker;

    // open file stream to look for points information
    std::ifstream fstreamp( filename.c_str() );
    if (fstreamp.fail()) {
        std::cerr << "Error in readNetgenMesh: File not found or locked" << std::endl;
        abort();
    }

    std::cout << "Reading netgen mesh file: " << filename << std::endl;
    getline(fstreamp, line);

    if( line.find("mesh3d") == std::string::npos ){
        std::cerr << "Error in readNetgenMesh: mesh file is not in mesh3d format (netgen)"\
                  << std::endl;
        abort();
    }

    /* I assume as I tested that faces stored are only boundary faces
       and edges stored are only a part of boundary ones with
       inside a face with a common marker, instead points can be
       inside the domain too, but this format doesn't say me
       which, so I'll find them from myself
    */
    flag=1|2|4|8;

    while( getline(fstreamp, line) ) {

        if( line.find("points") != std::string::npos && flag&1 ) {
            getline(fstreamp, line);
            std::stringstream parseline(line);

            parseline >> nVe;
            std::cout << "[readNetgenMesh] found " << nVe << " vertices... " << std::flush;

            bcnpoints.resize( nVe+1 );
            bpoints.resize( nVe+1 );
            pointcoor.resize( nDimensions * nVe );

            for(UInt i=0; i<nVe; i++) {
                // helping variables to read netgen file fields
                Real x, y, z;

                getline(fstreamp, line);
                std::stringstream parseline(line);

                parseline >> x >> y >> z;

                pointcoor[i*nDimensions] = x;
                pointcoor[i*nDimensions+1] = y;
                pointcoor[i*nDimensions+2] = z;
                bpoints[i] = false;
                bcnpoints[i] = NULLFLAG;
            }
            bpoints[nVe] = false;
            bcnpoints[nVe] = NULLFLAG;

            // done parsing point section
            flag&=~1;
            std::cout<< "loaded."<<std::endl;
            break;
        }
    }
    fstreamp.close();

    /* as Forma told, I have found sometimes problems with seekg, becouse
       it seems that sometimes it put the stream in a non consistent state,
       so I open new files instead of seeking, I have had not compile errors
       only sometimes problems on big mesh files with seekg
    */
    //  fstream.seekg(0,std::ios_base::beg);

    /* I assume as I tested that edges stored are only a part of
       boundary ones, so really they are meaningless to me
       unless I'm working with 2D meshes: in that case
       I can extract from here information about boundary vertices */
    std::ifstream fstreame(filename.c_str());
    while( getline(fstreame, line) ) {

        if( line.find("edgesegmentsgi2") != std::string::npos && flag&2 ) {
            getline(fstreame, line);
            std::stringstream parseline(line);

            parseline >> nBEd; // this will be discarded in 3D meshes

            std::cout << "[readNetgenMesh] found " << nBEd << " boundary edges... " << std::flush;
            edgepointID.resize(2*nBEd);

            for(UInt i=0; i<nBEd; ++i) {
                UInt surfnr, t, p1, p2;

                getline(fstreame, line);
                std::stringstream parseline(line);
                //          surfid  0   p1   p2   trignum1    trignum2   domin/surfnr1    domout/surfnr2
                //          ednr1   dist1   ednr2   dist2
                parseline >> surfnr >> t >> p1 >> p2; //>>t>>t>>t>>t>>t>>t>>t>>t;

                edgepointID[2*i] = p1;
                edgepointID[2*i+1] = p2; //>>t>>t>>t>>t>>t>>t>>t>>t;

            }

            // done parsing edge section
            flag&=~2;
            std::cout<< "loaded."<<std::endl;
            break;
        }
    }

    std::ifstream fstreamv(filename.c_str());
    while( getline(fstreamv, line) ) {

        if( line.find("volumeelements") != std::string::npos && flag&8) {
            getline(fstreamv, line);
            std::stringstream parseline(line);

            parseline >> nVo;
            std::cout << "[readNetgenMesh] found " << nVo << " volumes... " << std::flush;

            volumepointID.resize(4*nVo);

            for(UInt i=0; i<nVo; i++) {
                UInt matnr, np, p1, p2, p3, p4;

                getline(fstreamv, line);
                std::stringstream parseline(line);

                parseline >> matnr >> np >> p1 >> p2 >> p3 >> p4;

                volumepointID[4*i] = p1; volumepointID[4*i+1] = p2;
                volumepointID[4*i+2] = p3; volumepointID[4*i+3] = p4;

                ASSERT(np==4, "Error in readNetgenMesh: only tetrahedra elements supported")
                    }
            // done parsing volume section
            flag&=~8;
            std::cout<< "loaded."<<std::endl;
            break;
        }
    }
    fstreamv.close();

    std::ifstream fstreamf(filename.c_str());
    while( getline(fstreamf, line) ) {

        // TP 08/2008
        // surface elements section in a .vol file
        // is identified by key word "surfaceelements" in Netgen 4.5 (tested in RC2)
        // "surfaceelementsgi" is used in previous releases
        if( (line.find("surfaceelements") != std::string::npos) && flag&4) {
            getline(fstreamf, line);
            std::stringstream parseline(line);

            parseline >> nBFa;
            std::cout << "[readNetgenMesh] found " << nBFa
                      << " boundary faces... " << std::flush;

            facepointID.resize(3*nBFa);
            bcnsurf.resize(nBFa+1);
            bcnsurf[0]=NULLFLAG;

            for(UInt i=0; i<nBFa; i++) {
                UInt surfnr, bcnr, domin, domout, np, p1, p2, p3;

                getline(fstreamf, line);
                std::stringstream parseline(line);

                parseline >> surfnr >> bcnr >> domin >> domout >> np >> p1 >> p2 >> p3;

                facepointID[3*i] = p1;
                facepointID[3*i+1] = p2;
                facepointID[3*i+2] = p3;

                bcnsurf[i+1] = bcnr;
                //          std::cout<<"[readNetgenMesh] bcnr = " << bcnr << std::endl;
                ASSERT(np==3, "Error in readNetgenMesh: only triangular surfaces supported")

                //assume p1!=p2!=p3!=p1
                nBVe += bpoints[p1]?0:1;
                nBVe += bpoints[p2]?0:1;
                nBVe += bpoints[p3]?0:1;
                bpoints[p1]= bpoints[p2]= bpoints[p3]=true;

                /* here I set the boundary points marker
                   note: this works only with my patch
                   of strongerFlag
                   Face flag is assigned to face points
                   A point receives the "stronger" flag of the faces it belongs to
                    */
                bcnpoints[p1] = PMarker.setStrongerMarker(bcnpoints[p1],bcnr);
                bcnpoints[p2] = PMarker.setStrongerMarker(bcnpoints[p2],bcnr);
                bcnpoints[p3] = PMarker.setStrongerMarker(bcnpoints[p3],bcnr);

                /* now I have the surface and points, so I can calculate
                   the num of edges on the boundary, useful in case this is
                   a quadratic tetra, I calculate the bcnr too using
                   BareItemHandler<BareEdge> to store results
                */
                // the following is a complete list of edges in the 2D case
                // (only boundary edges in 3D)
                // may be useful even in the TWODIM case
                // (I've found this silly but easy way MM)
                BareEdge bed = setBareEdge(p1,p2);
                bihBedges.addIfNotThere(bed,(ID)NULLFLAG);
                bihBedges[bed]=(ID)EMarker.setStrongerMarker(bihBedges[bed], bcnr);

                bed=setBareEdge(p2,p3);
                bihBedges.addIfNotThere(bed,(ID)NULLFLAG);
                bihBedges[bed]=(ID)EMarker.setStrongerMarker(bihBedges[bed], bcnr);

                bed=setBareEdge(p3,p1);
                bihBedges.addIfNotThere(bed,(ID)NULLFLAG);
                bihBedges[bed]=(ID)EMarker.setStrongerMarker(bihBedges[bed], bcnr);

            }
            flag&=~4;
            // in the 3D case the only way to know the number of edges on boundary faces
            // is to count them!
            nBEd = bihBedges.howMany();
            std::cout<< "loaded."<<std::endl;
            break;
        }
    }
    fstreamf.close();

    ASSERT(flag==0, "[readNetgenMesh] the mesh file does not have all the required sections.")

    std::cout << "[readNetgenMesh] computed " << nBVe << " boundary vertices" << std::endl;

    // Euler formulas
    nFa=2*nVo+(nBFa/2);
    nEd=nVo+nVe+(3*nBFa-2*nBVe)/4;

    // Be a little verbose
    if ( GeoShape::numPoints > 4 ) {
        std::cout << "Quadratic Tetra  Mesh (from Linear geometry)" <<std::endl;
        nPo=nVe+nEd;
        nBPo=nBVe+nBEd;    // I calculated the real nBEd before
    } else {
        std::cout << "Linear Tetra Mesh" <<std::endl;
        nPo=nVe;
        nBPo=nBVe;
    }

    std::stringstream discardedLog;
    std::ostream& oStr = verbose ? std::cout : discardedLog;
    //points can be only vertices or on edges too

    std::cout << "#Vertices = "          << std::setw(10) << nVe
              << "  #BVertices       = " << std::setw(10) << nBVe << std::endl;
    oStr      << "#Faces    = "          << std::setw(10) << nFa << std::endl;
    oStr      << "  #Boundary Faces  = " << std::setw(10) << nBFa << std::endl;
    oStr      << "#Edges    = "          << std::setw(10) << nEd
              << "  #Boundary Edges  = " << std::setw(10) << nBEd << std::endl;
    std::cout << "#Points   = "          << std::setw(10) << nPo
              << "  #Boundary Points = " << std::setw(10) << nBPo << std::endl;
    std::cout << "#Volumes  = "          << std::setw(10) << nVo  << std::endl;

    // Set all basic data structure

    // I store all Points
    mesh.setMaxNumPoints   ( nPo, true );
    mesh.setMaxNumGlobalPoints( nPo );
    mesh.setNumBPoints     ( nBPo );
    mesh.setNumVertices    ( nVe );
    mesh.setNumGlobalVertices(nVe);
    mesh.setNumBVertices   ( nBVe );
    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges    ( nBEd );
    mesh.setMaxNumGlobalEdges ( nBEd );
    mesh.setNumEdges       ( nEd ); // Here the REAL number of edges (all of them)
    mesh.setNumBEdges      ( nBEd );
    // Only Boundary Faces
    mesh.setMaxNumFaces    ( nBFa );
    mesh.setMaxNumGlobalFaces ( nBFa );
    mesh.setNumFaces       ( nFa ); // Here the REAL number of faces (all of them)
    mesh.setNumBFaces      ( nBFa );

    mesh.setMaxNumVolumes  ( nVo, true );
    mesh.setMaxNumGlobalVolumes( nVo );

    mesh.setMarker         ( regionFlag ); // Add Marker to list of Markers

    typename RegionMesh3D<GeoShape,MC>::PointType * pp=0;
    typename RegionMesh3D<GeoShape,MC>::EdgeType * pe=0;
    typename RegionMesh3D<GeoShape,MC>::FaceType * pf=0;
    typename RegionMesh3D<GeoShape,MC>::VolumeType * pv=0;

    // addPoint()/Face()/Edge() returns a reference to the last stored point
    // I use that information to set all point info, by using a pointer.

    std::cout << "[readmesh3D] bpoints.size() = " << bpoints.size()
              << ", bcnpoints.size() = " << bcnpoints.size()
              << std::endl;

    for(UInt i=0; i<nVe; i++) {
        pp=&mesh.addPoint(bpoints[i+1]); //true if boundary point

        pp->setId     ( i + 1 );
        pp->setLocalId( i + 1 );
        mesh.localToGlobalNode().insert(std::make_pair(i+1, i+1));
        mesh.globalToLocalNode().insert(std::make_pair(i+1, i+1));

        pp->setMarker(bcnpoints[i+1]);
        pp->x()=pointcoor[nDimensions*i];
        pp->y()=pointcoor[nDimensions*i+1];
        pp->z()=pointcoor[nDimensions*i+2];
    }
    std::cout << "[readmesh3D] added points." << std::endl;

    /* here I set the real boundary edges that I stored
       in bihBedges
    */
    BareItemsHandler<BareEdge>::const_iterator bedge=bihBedges.begin();
    for(UInt i=0; i<nBEd; i++) {
    	UInt p1, p2;

        pe = &mesh.addEdge( true ); // Only boundary edges.
        pe->setMarker( EntityFlag(bedge->second) );
        p1=bedge->first.first;
        p2=bedge->first.second;
        pe->setPoint( 1, mesh.point( p1 ) ); // set edge conn.
        pe->setPoint( 2, mesh.point( p2 ) ); // set edge conn.
        bedge++;
    }
    std::cout << "[readmesh3D] added edges." << std::endl;

    for(UInt i=0;i<nVo;i++) {
    	UInt p1, p2, p3, p4;

        pv=&mesh.addVolume();
        pv->setId(i+1);
        pv->setLocalId(i+1);
        p1=volumepointID[4*i];
        p2=volumepointID[4*i+1];
        p3=volumepointID[4*i+2];
        p4=volumepointID[4*i+3];
        pv->setPoint(1, mesh.point(p1) );
        pv->setPoint(2, mesh.point(p2) );
        pv->setPoint(3, mesh.point(p3) );
        pv->setPoint(4, mesh.point(p4) );
    }
    std::cout << "[readmesh3D] added volumes." << std::endl;

    for(UInt i=0; i<nBFa; i++) {
    	UInt p1, p2, p3;

        pf=&mesh.addFace(true); // Only boundary faces
        p1=facepointID[3*i];
        p2=facepointID[3*i+1];
        p3=facepointID[3*i+2];

        pf->setMarker(EntityFlag(bcnsurf[i+1]));
        pf->setPoint(1,mesh.point(p1)); // set face conn.
        pf->setPoint(2,mesh.point(p2)); // set face conn.
        pf->setPoint(3,mesh.point(p3)); // set face conn.
    }
    std::cout << "[readmesh3D] added faces." << std::endl;

    // This part is to build a P2 mesh from a P1 geometry

    if ( GeoShape::numPoints > 4 ) p1top2(mesh);

    // Test mesh
    Switch sw;

    ///// CORRECTION JFG
    //if (mesh.check(1, true,true))done=0;

    if(!checkMesh3D(mesh, sw, true,verbose,oStr,std::cerr,oStr)) abort(); // CORRECTION JFG

    Real vols[3];
    getVolumeFromFaces(mesh, vols,oStr);
    oStr << "   VOLUME ENCLOSED BY THE MESH COMPUTED BY INTEGRATION ON"<<
        " BOUNDARY FACES"<<std::endl;
    oStr << "INT(X)     INT(Y)      INT(Z) <- they should be equal and equal to"
         << std::endl
         << "                                 the voulume enclosed by the mesh "
         << std::endl;
    oStr << vols[0] << " " << vols[1] << " " << vols[2] << std::endl;

    oStr << "   BOUNDARY FACES ARE DEFINING A CLOSED SURFACE IF "
         << testClosedDomain(mesh,oStr) << std::endl
         << " IS (ALMOST) ZERO" << std::endl;

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
