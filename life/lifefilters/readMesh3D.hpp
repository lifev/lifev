//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Mesh reader from mesh3d files

    @author Luca Formaggia <luca.formaggia@polimi.it>
    @contributor JFG, Nur Aiman Fadel <nur.fadel@mail.polimi.it>
    @maintainer Nur Aiman Fadel <nur.fadel@mail.polimi.it>

    @date 29-06-2002

    Mesh reader that it is able to read 3d meshes.<br>
    INRIAMesh used either spaces or CR as separators.<br>
 */

#ifndef _READMESH3D_HH_
#define _READMESH3D_HH_ 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>
#include <boost/lambda/lambda.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/util_string.hpp>

#include <life/lifemesh/bareItems.hpp>

#include <life/lifefilters/mesh_util.hpp>
#include <life/lifefilters/selectMarker.hpp>

namespace LifeV
{

struct FiveNumbers
{
public:
    UInt i1,i2,i3,i4;
    Int ibc;
};

// ===================================================
// Mpp mesh readers
// ===================================================

//! readMppFileHead - reads mesh++ Tetra meshes.
/*!
  It converts Tetra meshes into Quadratic Tetra if needed it.

  @param myStream,
  @param numberVertices,
  @param numberBoundaryVertices,
  @param numberBoundaryFaces,
  @param numberBoundaryEdges,
  @param numberVolumes,
  @return false if mesh check is unsuccessfull.
*/

bool
readMppFileHead( std::ifstream & myStream,
                 UInt          & numberVertices,
                 UInt          & numberBoundaryVertices,
                 UInt          & numberBoundaryFaces,
                 UInt          & numberBoundaryEdges,
                 UInt          & numberVolumes );

//! readMppFile - reads mesh++ Tetra meshes.
/*!
  It converts Tetra meshes into Quadratic Tetra if needed it.

  @param mesh, the mesh data structure to fill in.
  @param fileName, the name of the mesh file  to read.
  @param regionFlag, the identifier for the region.
  @param verbose, setting it as true, the output is verbose (the default is false).
  @return true if everything went fine, false otherwise.
*/

template <typename GeoShape, typename MC>
bool
readMppFile( RegionMesh3D<GeoShape, MC> & mesh,
             const std::string          & fileName,
             EntityFlag                   regionFlag,
             bool                         verbose = false )
{
    std::string line;

    Real x, y, z;

    Int ity, ity_id;

    UInt done = 0;
    UInt i;
    UInt nVe( 0 ), nBVe( 0 ), nFa( 0 ), nBFa( 0 ), nPo( 0 ), nBPo( 0 ), nEd( 0 ), nBEd( 0 );
    UInt nVo( 0 );
    UInt p1, p2, p3, p4;

    std::stringstream discardedLog;

    std::ostream& oStr = verbose ? std::cout : discardedLog;

    ASSERT_PRE0( GeoShape::Shape == TETRA ,  "readMppFiles reads only tetra meshes" ) ;
    ASSERT_PRE0( GeoShape::Shape == TETRA,   "Sorry, readMppFiles reads only tetra meshes" );
    ASSERT_PRE0( GeoShape::numVertices <= 6, "Sorry, readMppFiles handles only liner&quad tetras" );

    // open stream to read header

    std::ifstream hstream( fileName.c_str() );

    if ( hstream.fail() )
    {
        std::cerr << " Error in readMpp: File " << fileName
                  << " not found or locked" << std::endl;
        std::abort();
    }

    std::cout << "Reading mesh++ file" << std::endl;

    if ( ! readMppFileHead( hstream, nVe, nBVe, nBFa, nBEd, nVo ) )
    {
        std::cerr << " Error While reading mesh++ file headers" << std::endl;
        std::abort() ;
    }

    hstream.close();

    //Reopen the stream: I know it is stupid but this is how it goes
    std::ifstream myStream( fileName.c_str() );

    if ( myStream.fail() )
    {
        std::cerr << " Error in readMpp: File " << fileName
                  << " not found or locked" << std::endl;
        std::abort();
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
        nPo  = nVe;
        nBPo = nBVe;
    }

    std::cout << "Number of Vertices = "          << std::setw( 10 ) << nVe  << std::endl
              << "Number of Boundary Vertices = " << std::setw( 10 ) << nBVe << std::endl;
    oStr      << "Number of Faces    = "          << std::setw( 10 ) << nFa  << std::endl
              << "Number of Boundary Faces    = " << std::setw( 10 ) << nBFa << std::endl;
    oStr      << "Number of Edges    = "          << std::setw( 10 ) << nEd  << std::endl
              << "Number of Boundary Edges    = " << std::setw( 10 ) << nBEd << std::endl;
    std::cout << "Number of Points   = "          << std::setw( 10 ) << nPo  << std::endl
              << "Number of Boundary Points   = " << std::setw( 10 ) << nBPo << std::endl
              << "Number of Volumes  = "          << std::setw( 10 ) << nVo  << std::endl;

    // Set all basic data structure:

    // I store all Points
    mesh.setMaxNumPoints       ( nPo, true );
    mesh.setMaxNumGlobalPoints ( nPo );
    mesh.setNumBPoints         ( nBPo );
    mesh.setNumVertices        ( nVe );
    mesh.setNumGlobalVertices  (nVe);
    mesh.setNumBVertices       ( nBVe );

    // Only Boundary Edges
    mesh.setMaxNumEdges        ( nEd );
    mesh.setMaxNumGlobalEdges  ( nEd );
    mesh.setNumEdges           ( nEd ); // Here the REAL number of edges (all of them)
    mesh.setNumBEdges          ( nBEd );

    // Only Boundary Faces
    mesh.setMaxNumFaces        ( nBFa );
    mesh.setMaxNumGlobalFaces  ( nBFa );
    mesh.setNumFaces           ( nFa ); // Here the REAL number of edges (all of them)
    mesh.setNumBFaces          ( nBFa );

    mesh.setMaxNumVolumes      ( nVo, true );
    mesh.setMaxNumGlobalVolumes( nVo);

    mesh.setMarker             ( regionFlag ); // Mark the region


    typename RegionMesh3D<GeoShape, MC>::PointType  * pp = 0;
    typename RegionMesh3D<GeoShape, MC>::EdgeType   * pe = 0;
    typename RegionMesh3D<GeoShape, MC>::FaceType   * pf = 0;
    typename RegionMesh3D<GeoShape, MC>::VolumeType * pv = 0;

    // addPoint(), Face() and Edge() return a reference to the last stored point
    // I use that information to set all point info, by using a pointer.

    UInt count = 0;
    Int ibc;

    while ( next_good_line( myStream, line ).good() )
    {
        if ( line.find( "odes" ) != std::string::npos )
        {

            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );
            //      _numberVertices=atoi(node_s);

            for ( i = 0; i < nVe; i++ )
            {
#ifdef OLDMPPFILE
                myStream >> x >> y >> z >> ity >> ibc;
#else

                myStream >> x >> y >> z >> ity >> ity_id;
                if ( ity != 3 )
                    myStream >> ibc;
#endif

                if ( ity != 3 )
                {
                    ++count;
                    pp = &mesh.addPoint( true );
 
                  //Boundary point. Boundary switch set by the mesh method.

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

                mesh.localToGlobalNode().insert( std::make_pair( i + 1, i + 1) );
                mesh.globalToLocalNode().insert( std::make_pair( i + 1, i + 1) );
            }

            oStr << "Vertices Read " << std::endl;
            done++;

            if ( count != nBVe )
                std::cerr << "NumB points inconsistent!" << std::endl;
        }
        if ( line.find( "iangular" ) != std::string::npos )
        {
            oStr << "Reading Bfaces " << std::endl;
            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );

            for ( i = 0; i < nBFa; i++ )
            {

#ifdef OLDMPPFILE
                myStream >> p1 >> p2 >> p3 >> ity >> ibc;
#else

                myStream >> p1 >> p2 >> p3 >> ity >> ity_id >> ibc;
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
            //_numberBoundaryEdges=atoi(node_s);

            for ( i = 0; i < nBEd; i++ )
            {
#ifdef OLDMPPFILE
                myStream >> p1 >> p2 >> ity >> ibc;
#else
                myStream >> p1 >> p2 >> ity >> ity_id >> ibc;
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
                myStream >> p1 >> p2 >> p3 >> p4;
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
      {
        p1top2( mesh );
      }

    myStream.close();

    // Test mesh
    Switch sw;

    ///// CORRECTION JFG
    //if (mesh.check(1, true,true))done=0;

    if ( !checkMesh3D( mesh, sw, true, verbose, oStr, std::cerr, oStr ) )
      {
        std::abort(); // CORRECTION JFG
       }

    Real vols[ 3 ];
    getVolumeFromFaces( mesh, vols, oStr );

    oStr << "   VOLUME ENCLOSED BY THE MESH COMPUTED BY INTEGRATION ON"
         << " BOUNDARY FACES" << std::endl;
    oStr << "INT(X)     INT(Y)      INT(Z) <- they should be equal and equal to" << std::endl
         << "                                 the voulume enclosed by the mesh " << std::endl;
    oStr << vols[ 0 ] << " " << vols[ 1 ] << " " << vols[ 2 ] << std::endl;
    oStr << "   BOUNDARY FACES ARE DEFINING A CLOSED SURFACE IF "
         << testClosedDomain( mesh, oStr ) << std::endl
         << " IS (ALMOST) ZERO" << std::endl;

    return done == 4 ;
}// Function readMppFile

// ===================================================
// INRIA mesh readers
// ===================================================

//! nextIntINRIAMeshField -
/*!
  It gets an integer field from the std::string line
  if it is not empty, otherwise from the input stream.
  It assumes that the std::string is either empty or
  it contains and integer. No check is made to verify this.

  @param line, the mesh data structure to fill in.
  @param myStream, the name of the mesh file  to read.
  @return true if everything went fine, false otherwise.
*/

Int
nextIntINRIAMeshField( std::string const & line,
                       std::istream      & myStream );

//! readINRIAMeshFileHead - It Reads all basic info from INRIA MESH.
/*!
  It Reads all basic info from INRIA MESH file
  so as to be able to properly dimension all arrays

  @param myStream,
  @param numberVertices,
  @param numberBoundaryVertices,
  @param numberBoundaryFaces,
  @param numberBoundaryEdges,
  @param numberVolumes,
  @param numStoredFaces,
  @param shape,
  @param iSelect,
  @return false if mesh check is unsuccessfull.
*/

bool
readINRIAMeshFileHead( std::ifstream &        myStream,
                       UInt &                 numberVertices,
                       UInt &                 numberBoundaryVertices,
                       UInt &                 numberBoundaryFaces,
                       UInt &                 numberBoundaryEdges,
                       UInt &                 numberVolumes,
                       UInt &                 numStoredFaces,
                       ReferenceShapes &      shape,
                       InternalEntitySelector iSelect = InternalEntitySelector() );

//! readINRIAMeshFile - reads mesh++ Tetra meshes.
/*!
  It converts Tetra meshes into Quadratic Tetra if needed it.

  @param mesh, the mesh data structure to fill in.
  @param fileName, the name of the mesh file  to read.
  @param regionFlag, the identifier for the region.
  @param verbose, setting it as true, the output is verbose (the default is false).
  @param iSelect,
  @return true if everything went fine, false otherwise.
*/

template <typename GeoShape, typename MC>
bool
readINRIAMeshFile( RegionMesh3D<GeoShape, MC>&      mesh,
                   std::string const&               fileName,
                   EntityFlag                       regionFlag,
                   bool                             verbose = false,
                   InternalEntitySelector           iSelect = InternalEntitySelector() )
{
    std::string line;

    Real x, y, z;

    UInt done = 0;
    UInt i;
    UInt nVe( 0 ), nBVe( 0 ), nFa( 0 ), nBFa( 0 ), nPo( 0 ), nBPo( 0 ), nEd( 0 ), nBEd( 0 );
    UInt nVo( 0 );
    UInt numStoredFaces( 0 );
    UInt p1, p2, p3, p4, p5, p6, p7, p8;

    std::stringstream discardedLog;

    std::vector<FiveNumbers> faceHelp;

    ReferenceShapes shape( NONE );

    typename std::vector<FiveNumbers>::iterator faceHelpIterator;

    std::ostream& oStr = verbose ? std::cout : discardedLog;

    // open stream to read header

    std::ifstream hstream( fileName.c_str() );

    if ( verbose )
    {
        std::cout << "Reading form file " << fileName << std::endl;
    }

    if ( hstream.fail() )
    {
        std::cerr << " Error in readINRIAMeshFile: File " << fileName
                  << " not found or locked" << std::endl;
        std::abort();
    }

    if ( verbose )
    {
    std::cout << "Reading INRIA mesh file" << fileName << std::endl;
    }

    if ( ! readINRIAMeshFileHead( hstream, nVe, nBVe, nBFa, nBEd, nVo, numStoredFaces,shape,
                                  iSelect) )
    {
        std::cerr << " Error While reading INRIA mesh file headers" << std::endl;
        std::abort() ;
    }

    hstream.close();

    //Reopen the stream: I know it is stupid but this is how it goes
    std::ifstream myStream( fileName.c_str() );

    if ( myStream.fail() )
    {
        std::cerr << " Error in readINRIAMeshFile: File " << fileName
                  << " not found or locked" << std::endl;
        std::abort();
    }

    ASSERT_PRE0( GeoShape::Shape == shape, "INRIA Mesh file and mesh element shape is not consistent" );

    // Euler formulas to get number of faces and number of edges
    nFa = 2 * nVo + ( nBFa / 2 );
    Int num1  = nVe + nVo;
    Int num2  = nBVe;
    Int num3  = nBFa;

    nEd = ( 3 * num3 - 2 * num2 ) / 4 + num1;

//    nEd = (int) nVo + nVe + ( 3 * nBFa + dummy ) / 4;

    // Be a little verbose
    switch ( shape )
    {

    case HEXA:
        ASSERT_PRE0( GeoShape::numPoints == 8, "Sorry I can read only bilinear Hexa meshes" );
        std::cout << "Linear Hexa Mesh" << std::endl;
        nPo =  nVe;
        nBPo = nBVe;
        break;

    case TETRA:
        if ( GeoShape::numPoints > 4 )
        {
            //    if (GeoShape::numPoints ==6 )
            std::cout << "Quadratic Tetra  Mesh (from Linear geometry)" << std::endl;
            nPo = nVe + nEd;

            // nBPo=nBVe+nBEd; // FALSE : nBEd is not known at this stage in a INRIA file (JFG 07/2002)
            // I use the relation  nBVe + nBFa - 2 = nBEd, But, is it general (hole...)  (JFG 07/2002)

            nBEd = ( Int ( nBVe + nBFa - Int (2  ) ) > 0 ? ( nBVe + nBFa - 2 ) : 0 );
            nBPo = ( Int ( nBVe + ( nBVe + nBFa - Int ( 2 ) ) ) > 0 ?nBVe + ( nBVe + nBFa - 2 ) : 0 );
        }
        else
        {
            if ( verbose )
               {
                std::cout << "Linear Tetra Mesh" << std::endl;
               }

            nPo  = nVe;
            nBPo = nBVe;
            nBEd = ( Int ( nBVe + nBFa - Int( 2 ) ) > 0 ? ( nBVe + nBFa - 2 ) : 0 );
        }

        break;

    default:
        ERROR_MSG( "Current version of INRIA Mesh file reader only accepts TETRA and HEXA" );
    }

    oStr << "Number of Vertices        = "  << std::setw( 10 ) << nVe            << std::endl
         << "Number of BVertices       = "  << std::setw( 10 ) << nBVe           << std::endl
         << "Number of Faces           = "  << std::setw( 10 ) << nFa            << std::endl
         << "Number of Boundary Faces  = "  << std::setw( 10 ) << nBFa           << std::endl
         << "Number of Stored Faces    = "  << std::setw( 10 ) << numStoredFaces << std::endl
         << "Number of Edges           = "  << std::setw( 10 ) << nEd            << std::endl
         << "Number of Boundary Edges  = "  << std::setw( 10 ) << nBEd           << std::endl
         << "Number of Points          = "  << std::setw( 10 ) << nPo            << std::endl
         << "Number of Boundary Points = "  << std::setw( 10 ) << nBPo           << std::endl
         << "Number of Volumes         = "  << std::setw( 10 ) << nVo            << std::endl;

    // Set all basic data structure

    // I store all Points
    mesh.setMaxNumPoints       ( nPo, true );
    mesh.setMaxNumGlobalPoints ( nPo );
    mesh.setNumBPoints         ( nBPo );
    mesh.setNumVertices        ( nVe );
    mesh.setNumGlobalVertices  ( nVe );
    mesh.setNumBVertices       ( nBVe );
    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges        ( nBEd );
    mesh.setMaxNumGlobalEdges  ( nEd );
    mesh.setNumEdges           ( nEd ); // Here the REAL number of edges (all of them)
    mesh.setNumBEdges          ( nBEd );
    // Only Boundary Faces
    mesh.setMaxNumFaces        ( numStoredFaces );
    mesh.setMaxNumGlobalFaces  ( nBFa );
    mesh.setNumFaces           ( nFa ); // Here the REAL number of faces (all of them)
    mesh.setNumBFaces          ( nBFa );

    mesh.setMaxNumVolumes      ( nVo, true );
    mesh.setMaxNumGlobalVolumes( nVo );

    mesh.setMarker             ( regionFlag ); // Add Marker to list of Markers

    typedef typename RegionMesh3D<GeoShape, MC>::PointType  PointType;
    typedef typename RegionMesh3D<GeoShape, MC>::VolumeType VolumeType;


    typename RegionMesh3D<GeoShape, MC>::PointType  * pp = 0;
    typename RegionMesh3D<GeoShape, MC>::EdgeType   * pe = 0;
    typename RegionMesh3D<GeoShape, MC>::FaceType   * pf = 0;
    typename RegionMesh3D<GeoShape, MC>::VolumeType * pv = 0;
    // addPoint()/Face()/Edge() returns a reference to the last stored point
    // I use that information to set all point info, by using a pointer.

    UInt count = 0;
    Int  ibc;

    // To account for internal faces
    if ( numStoredFaces > nBFa )
    {
        faceHelp.resize( numStoredFaces - nBFa );
        faceHelpIterator = faceHelp.begin();

        oStr << "WARNING: The mesh file (apparently) contains "
             << numStoredFaces - nBFa << " internal faces" << std::endl;

    }

    while ( next_good_line( myStream, line ).good() )
    {
        if ( line.find( "Vertices" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), myStream );

            for ( i = 0; i < nVe; i++ )
            {
                myStream >> x >> y >> z >> ibc;

//                if (ibc == 1 ) ibc = 100;

                if ( !iSelect(EntityFlag(ibc)))
                {
                    ++count;
                // Boundary point. Boundary switch set by the mesh method.
                    pp = &mesh.addPoint( true ); 
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

                mesh.localToGlobalNode().insert( std::make_pair( i + 1, i + 1 ) );
                mesh.globalToLocalNode().insert( std::make_pair( i + 1, i + 1 ) );
            }

            oStr << "Vertices Read " << std::endl;
            oStr << "size of the node storage is "
                 << count * sizeof( PointType ) / 1024. / 1024. << std::endl;
            done++;

            if ( count != nBVe )
              {
                std::cerr << "NumB points inconsistent!" << std::endl;
              }
        }

        if ( line.find( "Triangles" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), myStream );
            oStr << "Reading Bfaces " << std::endl;

            for ( i = 0; i < numStoredFaces; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> ibc;

                if ( numStoredFaces > nBFa )
                {
                    if ( mesh.point( p1 ).boundary() && mesh.point( p2 ).boundary() &&
                            mesh.point( p3 ).boundary() )
                    {
                        pf = &( mesh.addFace( true ) ); // Boundary faces
                        pf->setMarker( EntityFlag( ibc ) );
                        pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                        pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                        pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.

                    }

                    else
                    {
                        faceHelpIterator->i1 = p1;
                        faceHelpIterator->i2 = p2;
                        faceHelpIterator->i3 = p3;
                        faceHelpIterator->ibc = ibc;

                        ++faceHelpIterator;
                    }
                }

                else
                {

                    pf = &( mesh.addFace( true ) ); // Only boundary faces

                    pf->setMarker( EntityFlag( ibc ) );
                    pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                    pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                    pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.
                }
            }

         for ( faceHelpIterator = faceHelp.begin(); faceHelpIterator != faceHelp.end(); ++faceHelpIterator )
          {
              p1  = faceHelpIterator->i1;
              p2  = faceHelpIterator->i2;
              p3  = faceHelpIterator->i3;
              ibc = faceHelpIterator->ibc;
              pf  = &( mesh.addFace( false ) ); // INTERNAL FACE
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
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), myStream );

            oStr << "Reading Bfaces " << std::endl;

            for ( i = 0; i < nBFa; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> p4 >> ibc;

                if ( numStoredFaces > nBFa )
                {
                    if ( mesh.point( p1 ).boundary() && mesh.point( p2 ).boundary() &&
                            mesh.point( p3 ).boundary() )
                    {
                        pf = &( mesh.addFace( true ) ); // Boundary faces
                        pf->setMarker( EntityFlag( ibc ) );
                        pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                        pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                        pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.
                        pf->setPoint( 4, mesh.point( p4 ) ); // set face conn.

                    }

                    else
                    {
                        faceHelpIterator->i1  = p1;
                        faceHelpIterator->i2  = p2;
                        faceHelpIterator->i3  = p3;
                        faceHelpIterator->i4  = p4;
                        faceHelpIterator->ibc = ibc;

                        ++faceHelpIterator;
                    }
                }

                else
                {
                    pf = &( mesh.addFace( true ) ); // Only boundary faces
                    pf->setMarker( EntityFlag( ibc ) );
                    pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                    pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                    pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.
                    pf->setPoint( 4, mesh.point( p4 ) ); // set face conn.
                }
            }

            oStr << "Boundary Faces Read " << std::endl;

            for ( faceHelpIterator=faceHelp.begin(); faceHelpIterator!=faceHelp.end(); ++faceHelpIterator )
            {
                p1 = faceHelpIterator->i1;
                p2 = faceHelpIterator->i2;
                p3 = faceHelpIterator->i3;
                p4 = faceHelpIterator->i4;
                ibc = faceHelpIterator->ibc;
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
            nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), myStream );
            oStr << "Reading Bedges " << std::endl;

            for ( i = 0; i < nBEd; i++ )
            {
                myStream >> p1 >> p2 >> ibc;
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
            nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), myStream );
            oStr << "Reading Volumes " << std::endl;

            for ( i = 0; i < nVo; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> p4 >> ibc;
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
            oStr << "size of the volume storage is " << sizeof( VolumeType ) * count / 1024. / 1024. 
                 << " Mo." << std::endl;
            oStr << count << " Volume elements Read" << std::endl;
            done++;
        }

        if ( line.find( "Hexahedra" ) != std::string::npos )
        {
            count = 0;
            nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), myStream );
            oStr << "Reading Volumes " << std::endl;
            for ( i = 0; i < nVo; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> p8 >> ibc;
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

    // Test the mesh
    Switch sw;

    // this if is the verbose version
    if ( !checkMesh3D( mesh, sw, true, verbose, oStr, std::cerr, oStr ) )
       {
        std::abort();
       }

    // This part is to build a P2 mesh from a P1 geometry
    if ( shape == TETRA && GeoShape::numPoints > 4 )
       {
        p1top2( mesh );
       }

    myStream.close();

    Real vols[ 3 ];
    getVolumeFromFaces( mesh, vols, oStr );

    oStr << "   VOLUME ENCLOSED BY THE MESH COMPUTED BY INTEGRATION ON" 
         << " BOUNDARY FACES" << std::endl;
    oStr << "INT(X)     INT(Y)      INT(Z) <- they should be equal and equal to" << std::endl 
         << "                                 the voulume enclosed by the mesh " << std::endl;
    oStr << vols[ 0 ] << " " << vols[ 1 ] << " " << vols[ 2 ] << std::endl;
    oStr << "   BOUNDARY FACES ARE DEFINING A CLOSED SURFACE IF " 
         << testClosedDomain( mesh, oStr ) << std::endl 
         << " IS (ALMOST) ZERO" << std::endl;

    return done == 4 ;
}// Function readINRIAMeshFile

// ===================================================
// GMSH mesh readers
// ===================================================

//! readGmshFile - it reads a GMSH mesh file
/*!
   It reads a 3D gmsh mesh file and store it in a RegionMesh3D.

   @param mesh mesh data structure to fill in
   @param fileName name of the gmsh mesh file  to read
   @param regionFlag identifier for the region
   @param verbose whether the function shall be verbose
   @return true if everything went fine, false otherwise
*/

template <typename GeoShape, typename MC>
bool
readGmshFile( RegionMesh3D<GeoShape, MC> & mesh,
              const std::string &          fileName,
              EntityFlag                   regionFlag,
              bool                         verbose = false )
{
    std::ifstream __is ( fileName.c_str() );

#ifdef DEBUG
    std::debug( 3000 ) << "gmsh reading: " << fileName << "\n";
#endif

    //    char __buf[256];
    std::string __buf;

    for (Int ii = 0; ii < 6; ++ii)
    {
        __is >> __buf;
        std::cout << "buf: " << __buf << "\n";
    }

    UInt __n;
    __is >> __n;

#ifdef DEBUG
    std::debug( 3000 ) << "number of nodes: " << __n;
#endif

    // Add Marker to list of Markers
    mesh.setMarker( regionFlag );


    std::vector<Real> __x( 3 * __n );
    std::vector<bool> __isonboundary( __n );
    std::vector<UInt> __whichboundary( __n );

#ifdef DEBUG
    std::debug( 3000 ) << "reading " << __n << " nodes\n";
#endif

    std::map<Int,Int> itoii;

    for ( UInt __i = 0; __i < __n; ++__i )
    {
        UInt __ni;
        __is >> __ni
        >> __x[ 3 * __i ]
        >> __x[ 3 * __i + 1 ]
        >> __x[ 3 * __i + 2 ];

        itoii[ __ni - 1] = __i;
    }
    __is >> __buf;

#ifdef DEBUG
    std::debug( 3000 ) << "buf: " << __buf << "\n";
#endif

    __is >> __buf;

#ifdef DEBUG
    std::debug( 3000 ) << "buf: " << __buf << "\n";
#endif

    UInt __nele;
    __is >> __nele;

    typename RegionMesh3D<GeoShape, MC>::EdgeType   * pe = 0;
    typename RegionMesh3D<GeoShape, MC>::FaceType   * pf = 0;
    typename RegionMesh3D<GeoShape, MC>::VolumeType * pv = 0;

#ifdef DEBUG
    std::debug( 3000 ) << "number of elements: " << __nele << "\n";
#endif

    std::vector<std::vector<int> > __e( __nele );
    std::vector<int>               __et( __nele );
    std::vector<int>               __etype( __nele );
    std::vector<int>               __gt( 32 );
    __gt.assign( 32, 0 );

    for ( UInt __i = 0; __i < __nele; ++__i )
    {
        Int __ne, __t, __np;

        //std::debug() << __i + 1 << " ";

        __is >> __buf;
        __is >> __ne;

        //std::debug() << __ne << " ";

        switch ( __ne )
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

#ifdef DEBUG
            std::debug( 3000 ) << "Element type unknown " << __ne << "\n";
#endif

            ASSERT( true, "Elements type unsupported.\n" )
        }

        __is >> __t;

        //std::debug() << __t << " ";

        bool ibcSet = false;
        Int  flag   = 0;
        Int __tag( 0 );

        for ( Int iflag = 0; iflag < __t; ++iflag )
        {
            __is >> flag;

            if ( !ibcSet )
            {
                __tag = flag;
                ibcSet = true;
            }
        }

        ++__gt[ __ne ];

        __etype[ __i ] = __ne;

        __et[ __i ] = __tag;
        __e[ __i ].resize( __np );

        Int __p = 0;

        while ( __p != __np )
        {
            Int node;
            __is >> node;
            __e[ __i ][ __p ] = node;
            __e[ __i ][ __p ] = itoii[ __e[ __i ][ __p ] - 1];
            __e[ __i ][ __p ] += 1;

            ++__p;
        }
    }


    // Euler formulas
    UInt n_volumes = __gt[ 4 ];
    UInt n_faces_boundary = __gt[ 2 ];
    UInt n_faces_total = 2 * n_volumes + ( n_faces_boundary / 2 );

    mesh.setMaxNumGlobalPoints( __n );
    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges( __gt[ 1 ] );
    mesh.setNumEdges   ( __gt[ 1 ] ); // Here the REAL number of edges (all of them)
    mesh.setNumBEdges  ( __gt[ 1 ] );
    mesh.setMaxNumGlobalEdges( __gt[ 1 ] );

#ifdef DEBUG
    std::debug( 3000 ) << "number of edges= " << __gt[ 1 ] << "\n";
#endif

    // Only Boundary Faces
    mesh.setMaxNumFaces( n_faces_total );
    mesh.setNumFaces   ( n_faces_total ); // Here the REAL number of edges (all of them)
    //mesh.setMaxNumFaces( n_faces_boundary );
    //mesh.setNumFaces   ( n_faces_boundary ); // Here the REAL number of edges (all of them)
    mesh.setNumBFaces  ( n_faces_boundary );
    mesh.setMaxNumGlobalFaces( n_faces_total );

#ifdef DEBUG
    std::debug( 3000 ) << "number of faces= " << n_faces_boundary << "\n";
#endif

    mesh.setMaxNumVolumes( n_volumes, true );
    mesh.setMaxNumGlobalVolumes( n_volumes );

#ifdef DEBUG
    std::debug( 3000 ) << "number of volumes= " << n_volumes << "\n";
#endif

    __isonboundary.assign( __n, false );
    __whichboundary.assign( __n, 0 );

    for ( UInt __i = 0; __i < __nele; ++__i )
    {
        switch ( __etype[ __i ] )
        {
            // triangular faces (linear)
        case 2:
        {
            __isonboundary[ __e[ __i ][ 0 ] - 1 ] = true;
            __isonboundary[ __e[ __i ][ 1 ] - 1 ] = true;
            __isonboundary[ __e[ __i ][ 2 ] - 1 ] = true;

            __whichboundary[ __e[ __i ][ 0 ] - 1 ] = __et[ __i ];
            __whichboundary[ __e[ __i ][ 1 ] - 1 ] = __et[ __i ];
            __whichboundary[ __e[ __i ][ 2 ] - 1 ] = __et[ __i ];
        }
        }
    }
    // add the point to the mesh
    typename RegionMesh3D<GeoShape, MC>::PointType * pp = 0;

    mesh.setMaxNumPoints( __n, true );
    mesh.setNumVertices ( __n );
    mesh.setNumBVertices( std::count( __isonboundary.begin(), __isonboundary.end(), true ) );
    mesh.setNumBPoints  ( mesh.numBVertices() );

#ifdef DEBUG
    std::debug( 3000 ) << "number of points : " << mesh.numPoints() << "\n";
    std::debug( 3000 ) << "number of boundary points : " << mesh.numBPoints() << "\n";
    std::debug( 3000 ) << "number of vertices : " << mesh.numVertices() << "\n";
    std::debug( 3000 ) << "number of boundary vertices : " << mesh.numBVertices() << "\n";
#endif

    for ( UInt __i = 0; __i < __n; ++__i )
    {
        pp = &mesh.addPoint( __isonboundary[ __i ] );
        pp->setMarker( __whichboundary[ __i ] );
        pp->setId     ( __i + 1 );
        pp->setLocalId( __i + 1 );
        pp->x() = __x[ 3 * __i ];
        pp->y() = __x[ 3 * __i + 1 ];
        pp->z() = __x[ 3 * __i + 2 ];
        mesh.localToGlobalNode().insert( std::make_pair( __i + 1, __i + 1 ) );
        mesh.globalToLocalNode().insert( std::make_pair( __i + 1, __i + 1 ) );
    }

    Int nVo = 1;
    // add the element to the mesh
    for ( UInt __i = 0; __i < __nele; ++__i )
    {
        switch ( __etype[ __i ] )
        {
        // segment(linear)
        case 1:
        {
            pe = &( mesh.addEdge( true ) );
            pe->setMarker( EntityFlag( __et[ __i ] ) );
            pe->setPoint( 1, mesh.point( __e[ __i ][ 0 ] ) );
            pe->setPoint( 2, mesh.point( __e[ __i ][ 1 ] ) );



        }
        break;

        // triangular faces (linear)
        case 2:
        {
            pf = &( mesh.addFace( true ) );
            pf->setMarker( EntityFlag( __et[ __i ] ) );
            pf->setPoint( 1, mesh.point( __e[ __i ][ 0 ] ) );
            pf->setPoint( 2, mesh.point( __e[ __i ][ 1 ] ) );
            pf->setPoint( 3, mesh.point( __e[ __i ][ 2 ] ) );

        }
        break;

        // quadrangular faces(linear)
        case 3:
        {
            pf = &( mesh.addFace( true ) );
            pf->setMarker( EntityFlag( __et[ __i ] ) );
            pf->setPoint( 1, mesh.point( __e[ __i ][ 0 ] ) );
            pf->setPoint( 2, mesh.point( __e[ __i ][ 1 ] ) );
            pf->setPoint( 3, mesh.point( __e[ __i ][ 2 ] ) );
            pf->setPoint( 4, mesh.point( __e[ __i ][ 3 ] ) );
        }
        break;

        // tetrahedrons(linear)
        case 4:
        {
            pv = &( mesh.addVolume() );
            pv->setId     ( nVo );
            pv->setLocalId( nVo++ );
//                pv->id() = __i + 1;
            pv->setMarker( EntityFlag( __et[ __i ] ) );
            pv->setPoint( 1, mesh.point( __e[ __i ][ 0 ] ) );
            pv->setPoint( 2, mesh.point( __e[ __i ][ 1 ] ) );
            pv->setPoint( 3, mesh.point( __e[ __i ][ 2 ] ) );
            pv->setPoint( 4, mesh.point( __e[ __i ][ 3 ] ) );
        }
        break;

        // hexahedrons(linear)
        case 5:
        {
            pv = &( mesh.addVolume() );

            pv->setId     ( __i + 1 );
            pv->setLocalId( __i + 1 );
//                pv->id() = __i + 1;
            pv->setMarker( EntityFlag( __et[ __i ] ) );
            pv->setPoint( 1, mesh.point( __e[ __i ][ 0 ] ) );
            pv->setPoint( 2, mesh.point( __e[ __i ][ 1 ] ) );
            pv->setPoint( 3, mesh.point( __e[ __i ][ 2 ] ) );
            pv->setPoint( 4, mesh.point( __e[ __i ][ 3 ] ) );
            pv->setPoint( 5, mesh.point( __e[ __i ][ 4 ] ) );
            pv->setPoint( 6, mesh.point( __e[ __i ][ 5 ] ) );
            pv->setPoint( 7, mesh.point( __e[ __i ][ 6 ] ) );
            pv->setPoint( 8, mesh.point( __e[ __i ][ 7 ] ) );
        }
        break;
        }
    }

    Switch sw;

    std::stringstream discardedLog;
    std::ostream& oStr = verbose ? std::cout : discardedLog;

    if ( checkMesh3D(mesh, sw, true,verbose,oStr,std::cerr,oStr) == false )
    {
        std::ostringstream __ex;
        __ex << "invalid mesh from GSMH";

        throw std::logic_error( __ex.str() );
    }

    return true;
} // Function readGmshFile

// ===================================================
// NetGen mesh readers
// ===================================================

//! readNetgenMesh - reads mesh++ Tetra meshes.
/*!
   It reads a 3D NetGen mesh file and store it in a RegionMesh3D.

   @param mesh mesh data structure to fill in.
   @param fileName name of the gmsh mesh file to read.
   @param regionFlag identifier for the region.
   @param verbose whether the function shall be verbose.
   @return true if everything went fine, false otherwise.
*/

template<typename GeoShape, typename MC>
bool
readNetgenMesh(RegionMesh3D<GeoShape,MC> & mesh,
               const std::string  &        fileName,
               EntityFlag                  regionFlag,
               bool                        verbose = false )
{
    // I will extract lines from iostream
    std::string line;

    // number of Geo Elements
    UInt nVe( 0 ), nBVe( 0 ), nPo( 0 ), nBPo( 0 );
    UInt nEd( 0 ), nBEd( 0 ), nFa( 0 ), nBFa( 0 );
    UInt nVo( 0 );

    // During the first access to file, build list of structures
    Vector pointcoor;
    std::vector<UInt> edgepointID, facepointID, volumepointID;

    // for poit i-th: is it a boundary point?
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
    std::ifstream fstreamp( fileName.c_str() );
    if (fstreamp.fail())
    {
        std::cerr << "Error in readNetgenMesh: File not found or locked" << std::endl;
        std::abort();
    }

    std::cout << "Reading netgen mesh file: " << fileName << std::endl;
    getline( fstreamp, line );

    if ( line.find( "mesh3d" ) == std::string::npos )
    {
        std::cerr << "Error in readNetgenMesh: mesh file is not in mesh3d format (netgen)"\
                  << std::endl;
        std::abort();
    }

    /*
       I assume as I tested that faces stored are only boundary faces
       and edges stored are only a part of boundary ones with
       inside a face with a common marker, instead points can be
       inside the domain too, but this format doesn't say me
       which, so I'll find them from myself
    */

    flag=1|2|4|8;

    while ( getline( fstreamp, line ) )
    {

        if ( line.find( "points" ) != std::string::npos && flag&1 )
        {
            getline( fstreamp, line );
            std::stringstream parseline( line );

            parseline >> nVe;
            std::cout << "[readNetgenMesh] found " << nVe << " vertices... " << std::flush;

            bcnpoints.resize( nVe + 1 );
            bpoints.resize( nVe + 1 );
            pointcoor.resize( nDimensions * nVe );

            for ( UInt i = 0; i < nVe; i++ )
            {
                // helping variables to read netgen file fields
                Real x, y, z;

                getline( fstreamp, line );
                std::stringstream parseline( line );

                parseline >> x >> y >> z;

                pointcoor[ i * nDimensions ]     = x;
                pointcoor[ i * nDimensions + 1 ] = y;
                pointcoor[ i * nDimensions + 2 ] = z;
                bpoints  [ i ] = false;
                bcnpoints[ i ] = NULLFLAG;
            }
            bpoints  [ nVe ] = false;
            bcnpoints[ nVe ] = NULLFLAG;

            // done parsing point section
            flag&=~1;
            std::cout << "loaded." << std::endl;
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
    std::ifstream fstreame( fileName.c_str() );
    while ( getline( fstreame, line ) )
    {

        if ( line.find( "edgesegmentsgi2" ) != std::string::npos && flag&2 )
        {
            getline( fstreame, line );
            std::stringstream parseline( line );

            parseline >> nBEd; // this will be discarded in 3D meshes

            std::cout << "[readNetgenMesh] found " << nBEd << " boundary edges... " << std::flush;

            edgepointID.resize( 2 * nBEd );

            for ( UInt i = 0; i < nBEd; ++i)
            {
                UInt surfnr, t, p1, p2;

                getline( fstreame, line );
                std::stringstream parseline( line );

                parseline >> surfnr >> t >> p1 >> p2;

                edgepointID[ 2 * i ] = p1;
                edgepointID[ 2 * i + 1 ] = p2;

            }

            // done parsing edge section
            flag&=~2;
            std::cout << "loaded." << std::endl;
            break;
        }
    }

    std::ifstream fstreamv( fileName.c_str() );
    while ( getline( fstreamv, line ) )
    {

        if ( line.find( "volumeelements" ) != std::string::npos && flag&8 )
        {
            getline( fstreamv, line );
            std::stringstream parseline( line );

            parseline >> nVo;
            std::cout << "[readNetgenMesh] found " << nVo << " volumes... " << std::flush;

            volumepointID.resize( 4 * nVo );

            for ( UInt i = 0; i < nVo; i++ )
            {
                UInt matnr, np, p1, p2, p3, p4;

                getline( fstreamv, line );
                std::stringstream parseline( line );

                parseline >> matnr >> np >> p1 >> p2 >> p3 >> p4;

                volumepointID[ 4 * i ] = p1;
                volumepointID[ 4 * i + 1 ] = p2;
                volumepointID[ 4 * i + 2 ] = p3;
                volumepointID[ 4 * i + 3 ] = p4;

                ASSERT( np==4, "Error in readNetgenMesh: only tetrahedra elements supported" )
            }
            // done parsing volume section
            flag&=~8;
            std::cout << "loaded." << std::endl;
            break;
        }
    }
    fstreamv.close();

    std::ifstream fstreamf( fileName.c_str() );
    while ( getline( fstreamf, line ) )
    {

        // TP 08/2008
        // surface elements section in a .vol file
        // is identified by key word "surfaceelements" in Netgen 4.5 (tested in RC2)
        // "surfaceelementsgi" is used in previous releases
        if ( ( line.find("surfaceelements") != std::string::npos ) && flag&4 )
        {
            getline( fstreamf, line );
            std::stringstream parseline( line );

            parseline >> nBFa;

            std::cout << "[readNetgenMesh] found " << nBFa
                      << " boundary faces... "     << std::flush;

            facepointID.resize( 3 * nBFa );
            bcnsurf.resize( nBFa + 1 );
            bcnsurf[ 0 ] = NULLFLAG;

            for ( UInt i = 0; i < nBFa; i++ )
            {
                UInt surfnr, bcnr, domin, domout, np, p1, p2, p3;

                getline( fstreamf, line );

                if ( line.empty() )
                {
                    // newer version of netgen inserts a blank line
                    // after each line in "surface elements" section
                    // make sure we ignore it
//                    std::cout << "\nfound empty line in netgen file!" << std::endl;
                    getline( fstreamf, line );
                }
                std::stringstream parseline( line );

                parseline >> surfnr >> bcnr >> domin >> domout >> np >> p1 >> p2 >> p3;

                facepointID[ 3 * i ] = p1;
                facepointID[ 3 * i + 1 ] = p2;
                facepointID[ 3 * i + 2 ] = p3;

                bcnsurf[ i + 1 ] = bcnr;
                //          std::cout<<"[readNetgenMesh] bcnr = " << bcnr << std::endl;
                ASSERT( np==3, "Error in readNetgenMesh: only triangular surfaces supported" )

                //assume p1!=p2!=p3!=p1
                nBVe += bpoints[ p1 ] ? 0:1;
                nBVe += bpoints[ p2 ] ? 0:1;
                nBVe += bpoints[ p3 ] ? 0:1;
                bpoints[ p1 ] = bpoints[ p2 ] = bpoints[ p3 ] = true;

                /* here I set the boundary points marker
                   note: this works only with my patch
                   of strongerFlag
                   Face flag is assigned to face points
                   A point receives the "stronger" flag of the faces it belongs to
                    */
                bcnpoints[ p1 ] = PMarker.setStrongerMarker( bcnpoints[ p1 ], bcnr );
                bcnpoints[ p2 ] = PMarker.setStrongerMarker( bcnpoints[ p2 ], bcnr );
                bcnpoints[ p3 ] = PMarker.setStrongerMarker( bcnpoints[ p3 ], bcnr );

                /* now I have the surface and points, so I can calculate
                   the num of edges on the boundary, useful in case this is
                   a quadratic tetra, I calculate the bcnr too using
                   BareItemHandler<BareEdge> to store results
                */
                // the following is a complete list of edges in the 2D case
                // (only boundary edges in 3D)
                // may be useful even in the TWODIM case
                // (I've found this silly but easy way MM)
                BareEdge bed = setBareEdge( p1, p2 );
                bihBedges.addIfNotThere( bed, ( ID )NULLFLAG );
                bihBedges[ bed ] = ( ID )EMarker.setStrongerMarker( bihBedges[ bed ], bcnr );

                bed = setBareEdge( p2, p3 );
                bihBedges.addIfNotThere( bed,( ID )NULLFLAG );
                bihBedges[ bed ] = ( ID )EMarker.setStrongerMarker( bihBedges[ bed ], bcnr );

                bed = setBareEdge( p3, p1 );
                bihBedges.addIfNotThere( bed,( ID )NULLFLAG );
                bihBedges[ bed ] = ( ID )EMarker.setStrongerMarker( bihBedges[ bed ], bcnr );

            }
            flag&=~4;
            // in the 3D case the only way to know the number of edges on boundary faces
            // is to count them!
            nBEd = bihBedges.howMany();
            std::cout << "loaded." << std::endl;
            break;
        }
    }

    fstreamf.close();

    ASSERT( flag==0, "[readNetgenMesh] the mesh file does not have all the required sections." )

    std::cout << "[readNetgenMesh] computed " << nBVe << " boundary vertices" << std::endl;

    // Euler formulas
    nFa = 2 * nVo + ( nBFa / 2 );
    nEd = nVo + nVe + ( 3 * nBFa - 2 * nBVe ) / 4;

    // Be a little verbose
    if ( GeoShape::numPoints > 4 )
    {
        std::cout << "Quadratic Tetra  Mesh (from Linear geometry)" <<std::endl;
        nPo = nVe + nEd;
        nBPo = nBVe + nBEd;    // I calculated the real nBEd before
    }

    else
    {
        std::cout << "Linear Tetra Mesh" <<std::endl;
        nPo = nVe;
        nBPo = nBVe;
    }

    std::stringstream discardedLog;
    std::ostream& oStr = verbose ? std::cout : discardedLog;
    //points can be only vertices or on edges too

    std::cout << "Number of Vertices        = " << std::setw( 10 ) << nVe  << std::endl
              << "Number of BVertices       = " << std::setw( 10 ) << nBVe << std::endl;
    oStr      << "Number of Faces           = " << std::setw( 10 ) << nFa  << std::endl
              << "Number of Boundary Faces  = " << std::setw( 10 ) << nBFa << std::endl
              << "Number of Edges           = " << std::setw( 10 ) << nEd  << std::endl
              << "Number of Boundary Edges  = " << std::setw( 10 ) << nBEd << std::endl;
    std::cout << "Number of Points          = " << std::setw( 10 ) << nPo  << std::endl
              << "Number of Boundary Points = " << std::setw( 10 ) << nBPo << std::endl
              << "Number of Volumes         = " << std::setw( 10 ) << nVo  << std::endl;

    // Set all basic data structure

    // I store all Points
    mesh.setMaxNumPoints       ( nPo, true );
    mesh.setMaxNumGlobalPoints ( nPo );
    mesh.setNumBPoints         ( nBPo );
    mesh.setNumVertices        ( nVe );
    mesh.setNumGlobalVertices  ( nVe );
    mesh.setNumBVertices       ( nBVe );

    // Only Boundary Edges
    mesh.setMaxNumEdges        ( nBEd );
    mesh.setMaxNumGlobalEdges  ( nBEd );
    mesh.setNumEdges           ( nEd ); // Here the REAL number of edges (all of them)
    mesh.setNumBEdges          ( nBEd );

    // Only Boundary Faces
    mesh.setMaxNumFaces        ( nBFa );
    mesh.setMaxNumGlobalFaces  ( nBFa );
    mesh.setNumFaces           ( nFa ); // Here the REAL number of faces (all of them)
    mesh.setNumBFaces          ( nBFa );

    mesh.setMaxNumVolumes      ( nVo, true );
    mesh.setMaxNumGlobalVolumes( nVo );

    mesh.setMarker             ( regionFlag ); // Add Marker to list of Markers

    typename RegionMesh3D<GeoShape,MC>::PointType  * pp = 0;
    typename RegionMesh3D<GeoShape,MC>::EdgeType   * pe = 0;
    typename RegionMesh3D<GeoShape,MC>::FaceType   * pf = 0;
    typename RegionMesh3D<GeoShape,MC>::VolumeType * pv = 0;

    // addPoint()/Face()/Edge() returns a reference to the last stored point
    // I use that information to set all point info, by using a pointer.

    std::cout << "[readmesh3D] bpoints.size() = " << bpoints.size()
              << ", bcnpoints.size() = "          << bcnpoints.size()
              << std::endl;

    for ( UInt i = 0; i < nVe; i++ )
    {
        pp=&mesh.addPoint( bpoints[ i + 1 ] ); //true if boundary point

        pp->setId     ( i + 1 );
        pp->setLocalId( i + 1 );
        mesh.localToGlobalNode().insert( std::make_pair( i + 1, i + 1 ) );
        mesh.globalToLocalNode().insert( std::make_pair( i + 1, i + 1 ) );

        pp->setMarker( bcnpoints[ i + 1 ] );
        pp->x() = pointcoor[ nDimensions * i ];
        pp->y() = pointcoor[ nDimensions * i + 1 ];
        pp->z() = pointcoor[ nDimensions * i + 2 ];
    }

    std::cout << "[readmesh3D] added points." << std::endl;

    /*
       here I set the real boundary edges that I stored
       in bihBedges
    */
    
    BareItemsHandler<BareEdge>::const_iterator bedge = bihBedges.begin();

    for ( UInt i=0; i < nBEd; i++ )
    {
        UInt p1, p2;

        pe = &mesh.addEdge( true ); // Only boundary edges.
        pe->setMarker( EntityFlag( bedge->second ) );
        p1 = bedge->first.first;
        p2 = bedge->first.second;
        pe->setPoint( 1, mesh.point( p1 ) ); // set edge conn.
        pe->setPoint( 2, mesh.point( p2 ) ); // set edge conn.

        bedge++;
    }

    std::cout << "[readmesh3D] added edges." << std::endl;

    for ( UInt i = 0; i < nVo; i++ )
    {
        UInt p1, p2, p3, p4;

        pv = &mesh.addVolume();
        pv->setId( i + 1 );
        pv->setLocalId(i + 1 );
        p1 = volumepointID[ 4 * i ];
        p2 = volumepointID[ 4 * i + 1 ];
        p3 = volumepointID[ 4 * i + 2 ];
        p4 = volumepointID[ 4 * i + 3 ];
        pv->setPoint( 1, mesh.point( p1 ) );
        pv->setPoint( 2, mesh.point( p2 ) );
        pv->setPoint( 3, mesh.point( p3 ) );
        pv->setPoint( 4, mesh.point( p4 ) );
    }

    std::cout << "[readmesh3D] added volumes." << std::endl;

    for ( UInt i = 0; i < nBFa; i++ )
    {
        UInt p1, p2, p3;

        pf = &mesh.addFace( true ); // Only boundary faces
        p1 = facepointID[ 3 * i ];
        p2 = facepointID[ 3 * i + 1 ];
        p3 = facepointID[ 3 * i + 2 ];

        pf->setMarker( EntityFlag( bcnsurf[ i + 1 ] ) );
        pf->setPoint( 1, mesh.point( p1 ) ); // set face conn.
        pf->setPoint( 2, mesh.point( p2 ) ); // set face conn.
        pf->setPoint( 3, mesh.point( p3 ) ); // set face conn.
    }
    std::cout << "[readmesh3D] added faces." << std::endl;

    // This part is to build a P2 mesh from a P1 geometry

    if ( GeoShape::numPoints > 4 ) 
    {
    p1top2( mesh );
    }

    // Test mesh
    Switch sw;

    ///// CORRECTION JFG
    //if (mesh.check(1, true,true))done=0;

    if ( !checkMesh3D( mesh, sw, true, verbose, oStr, std::cerr, oStr ) ) 
    {
    std::abort(); // CORRECTION JFG
    }

    Real vols[ 3 ];

    getVolumeFromFaces( mesh, vols,oStr );

    oStr << "   VOLUME ENCLOSED BY THE MESH COMPUTED BY INTEGRATION ON" << " BOUNDARY FACES" << std::endl
         << "INT(X)     INT(Y)      INT(Z) <- they should be equal and equal to" << std::endl
         << "                                 the voulume enclosed by the mesh " << std::endl
         << vols[ 0 ] << " " << vols[ 1 ] << " " << vols[ 2 ] << std::endl
         << "   BOUNDARY FACES ARE DEFINING A CLOSED SURFACE IF " 
         << testClosedDomain( mesh, oStr ) << std::endl
         << " IS (ALMOST) ZERO" << std::endl;

    return true;
} // Function readNetgenMesh

//! saveNetgenSolution -
/*!
  Ripped "from src/ng431/libsrc/interface/importsolution.cpp"

  @param fileName, the name of the mesh file  to read.
  @param U,
  @param fctname, default to "u".
  @return void.
*/

template <typename VectorType>
void
saveNetgenSolution(std::string       fileName,
                   const VectorType& solution,
                   std::string       functionName = "u")
{
    std::ofstream of(fileName.c_str());

    of << "solution " << functionName << " -size =" << solution.size() << " -components=1 -type=nodal" << std::endl;

    for ( UInt i = 0; i < solution.size(); ++i )
    {
        of << solution( i ) << std::endl;
    }

    of.close();
} // Function saveNetgenSolution

} // Namespace LifeV

#endif /* READMESH3D_H */
