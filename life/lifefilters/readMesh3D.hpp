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
#include <boost/numeric/ublas/vector.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/util_string.hpp>

#include <life/lifemesh/bareItems.hpp>

#include <life/lifefilters/mesh_util.hpp>
#include <life/lifefilters/selectMarker.hpp>

namespace LifeV
{

// @name Public typedefs
//@{
typedef boost::numeric::ublas::vector<Real> Vector;
//@}

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
             entityFlag_Type              regionFlag,
             bool                         verbose = false )
{
    std::string line;

    Real x, y, z;

    Int ity, ity_id;

    UInt done = 0;
    UInt i;
    UInt numberVertices( 0 ), numberBoundaryVertices( 0 ),
         numberFaces   ( 0 ), numberBoundaryFaces( 0 ),
         numberPoints  ( 0 ), numberBoundaryPoints( 0 ),
         numberEdges   ( 0 ), numberBoundaryEdges( 0 );
    UInt numberVolumes ( 0 );
    UInt p1, p2, p3, p4;

    std::stringstream discardedLog;

    std::ostream& oStr = verbose ? std::cout : discardedLog;

    ASSERT_PRE0( GeoShape::S_shape == TETRA ,  "readMppFiles reads only tetra meshes" ) ;
    ASSERT_PRE0( GeoShape::S_shape == TETRA,   "Sorry, readMppFiles reads only tetra meshes" );
    ASSERT_PRE0( GeoShape::S_numVertices <= 6, "Sorry, readMppFiles handles only liner&quad tetras" );

    // open stream to read header

    std::ifstream hstream( fileName.c_str() );

    if ( hstream.fail() )
    {
        std::cerr << " Error in readMpp: File " << fileName
                  << " not found or locked" << std::endl;
        std::abort();
    }

    std::cout << "Reading mesh++ file" << std::endl;

    if ( ! readMppFileHead( hstream, numberVertices, numberBoundaryVertices,
                            numberBoundaryFaces, numberBoundaryEdges, numberVolumes ) )
    {
        std::cerr << " Error While reading mesh++ file headers" << std::endl;
        std::abort() ;
    }

    hstream.close();

    //Reopen the stream: I know it is stupid but this is how it goes
    std::ifstream myStream( fileName.c_str() );

    if ( myStream.fail() )
    {
        std::cerr << " Error in readMpp: file " << fileName
                  << " not found or locked" << std::endl;
        std::abort();
    }

    // Euler formulas
    numberFaces = 2 * numberVolumes + ( numberBoundaryFaces / 2 );
    numberEdges = numberVolumes + numberVertices +
                  ( 3 * numberBoundaryFaces - 2 * numberBoundaryVertices ) / 4;

    // Be a little verbose
    if ( GeoShape::S_numPoints > 4 )
    {
        std::cout << "Quadratic Tetra  Mesh (from Linear geometry)"
                  << std::endl;
        numberPoints         = numberVertices         + numberEdges;
        numberBoundaryPoints = numberBoundaryVertices + numberBoundaryEdges;
    }

    else
    {
        std::cout << "Linear Tetra Mesh" << std::endl;
        numberPoints         = numberVertices;
        numberBoundaryPoints = numberBoundaryVertices;
    }

    std::cout << "Number of Vertices          = " << std::setw( 10 ) << numberVertices  << std::endl
              << "Number of Boundary Vertices = " << std::setw( 10 ) << numberBoundaryVertices << std::endl;
    oStr      << "Number of Faces             = " << std::setw( 10 ) << numberFaces  << std::endl
              << "Number of Boundary Faces    = " << std::setw( 10 ) << numberBoundaryFaces << std::endl
              << "Number of Edges             = " << std::setw( 10 ) << numberEdges  << std::endl
              << "Number of Boundary Edges    = " << std::setw( 10 ) << numberBoundaryEdges << std::endl;
    std::cout << "Number of Points            = " << std::setw( 10 ) << numberPoints << std::endl
              << "Number of Boundary Points   = " << std::setw( 10 ) << numberBoundaryPoints << std::endl
              << "Number of Volumes           = " << std::setw( 10 ) << numberVolumes  << std::endl;

    // Set all basic data structure:

    // I store all Points
    mesh.setMaxNumPoints       ( numberPoints, true );
    mesh.setMaxNumGlobalPoints ( numberPoints );
    mesh.setNumBPoints         ( numberBoundaryPoints );
    mesh.setNumVertices        ( numberVertices );
    mesh.setNumGlobalVertices  ( numberVertices );
    mesh.setNumBVertices       ( numberBoundaryVertices );

    // Only Boundary Edges
    mesh.setMaxNumEdges        ( numberEdges );
    mesh.setMaxNumGlobalEdges  ( numberEdges );
    mesh.setNumEdges           ( numberEdges ); // Here the REAL number of edges (all of them)
    mesh.setNumBEdges          ( numberBoundaryEdges );

    // Only Boundary Faces
    mesh.setMaxNumFaces        ( numberBoundaryFaces );
    mesh.setMaxNumGlobalFaces  ( numberBoundaryFaces );
    mesh.setNumFaces           ( numberFaces ); // Here the REAL number of edges (all of them)
    mesh.setNumBFaces          ( numberBoundaryFaces );

    mesh.setMaxNumVolumes      ( numberVolumes, true );
    mesh.setMaxNumGlobalVolumes( numberVolumes );

    mesh.setMarker             ( regionFlag ); // Mark the region


    typename RegionMesh3D<GeoShape, MC>::PointType  * pointerPoint  = 0;
    typename RegionMesh3D<GeoShape, MC>::EdgeType   * pointerEdge   = 0;
    typename RegionMesh3D<GeoShape, MC>::FaceType   * pointerFace   = 0;
    typename RegionMesh3D<GeoShape, MC>::VolumeType * pointerVolume = 0;

    // addPoint(), Face() and Edge() return a reference to the last stored point
    // I use that information to set all point info, by using a pointer.

    UInt count = 0;
    Int ibc;

    while ( nextGoodLine( myStream, line ).good() )
    {
        if ( line.find( "odes" ) != std::string::npos )
        {

            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );
            //      _numberVertices=atoi(node_s);

            for ( i = 0; i < numberVertices; i++ )
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
                    pointerPoint = &mesh.addPoint( true );

                  //Boundary point. Boundary switch set by the mesh method.

                    pointerPoint->setMarker( entityFlag_Type( ibc ) );
                }
                else
                {
                    pointerPoint = &mesh.addPoint( false );
                }

                pointerPoint->x() = x;
                pointerPoint->y() = y;
                pointerPoint->z() = z;
                pointerPoint->setMarker( entityFlag_Type( ibc ) );

                pointerPoint->setId     ( i + 1 );
                pointerPoint->setLocalId( i + 1 );

                mesh.localToGlobalNode().insert( std::make_pair( i + 1, i + 1) );
                mesh.globalToLocalNode().insert( std::make_pair( i + 1, i + 1) );
            }

            oStr << "Vertices Read " << std::endl;
            done++;

            if ( count != numberBoundaryVertices )
                std::cerr << "NumB points inconsistent!" << std::endl;
        }
        if ( line.find( "iangular" ) != std::string::npos )
        {
            oStr << "Reading Bfaces " << std::endl;
            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );

            for ( i = 0; i < numberBoundaryFaces; i++ )
            {

#ifdef OLDMPPFILE
                myStream >> p1 >> p2 >> p3 >> ity >> ibc;
#else

                myStream >> p1 >> p2 >> p3 >> ity >> ity_id >> ibc;
#endif

                pointerFace = &( mesh.addFace( true ) ); // Only boundary faces

                pointerFace->setMarker( entityFlag_Type( ibc ) );
                pointerFace->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                pointerFace->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                pointerFace->setPoint( 3, mesh.point( p3 ) ); // set face conn.
            }

            oStr << "Boundary faces read " << std::endl;
            done++;
        }

        if ( line.find( "Sides" ) != std::string::npos )
        {
            oStr << "Reading boundary edges " << std::endl;
            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );
            //_numberBoundaryEdges=atoi(node_s);

            for ( i = 0; i < numberBoundaryEdges; i++ )
            {
#ifdef OLDMPPFILE
                myStream >> p1 >> p2 >> ity >> ibc;
#else
                myStream >> p1 >> p2 >> ity >> ity_id >> ibc;
#endif
                pointerEdge = &mesh.addEdge( true ); // Only boundary edges.
                pointerEdge->setMarker( entityFlag_Type( ibc ) );
                pointerEdge->setPoint( 1, mesh.point( p1 ) ); // set edge conn.
                pointerEdge->setPoint( 2, mesh.point( p2 ) ); // set edge conn.
            }

            oStr << "Boundary edges read " << std::endl;
            done++;
        }

        count = 0;

        if ( line.find( "etrahedral" ) != std::string::npos )
        {
            oStr << "Reading volumes " << std::endl;
            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );

            for ( i = 0; i < numberVolumes; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> p4;
                pointerVolume = &mesh.addVolume();
                pointerVolume->setId     ( i + 1 );
                pointerVolume->setLocalId( i + 1);
//                pointerVolume->id() = i + 1;
                pointerVolume->setPoint( 1, mesh.point( p1 ) );
                pointerVolume->setPoint( 2, mesh.point( p2 ) );
                pointerVolume->setPoint( 3, mesh.point( p3 ) );
                pointerVolume->setPoint( 4, mesh.point( p4 ) );
                count++;
            }

            oStr << count << " Volume elements read" << std::endl;
            done++;
        }
    }
    // This part is to build a P2 mesh from a P1 geometry

    if ( GeoShape::S_numPoints > 4 )
    {
    	MeshUtility::p2MeshFromP1Data( mesh );
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

    oStr << "Volume enclosed by the mesh computed by integration on boundary faces "  << std::endl;
    oStr << "INT(X)     INT(Y)      INT(Z)      <- they should be equal and equal to" << std::endl
         << "                                   the voulume enclosed by the mesh "    << std::endl;
    oStr << vols[ 0 ] << "      " << vols[ 1 ] << "      " << vols[ 2 ] << std::endl;
    oStr << "Boundary faces are defining a closed surface if "
         << testClosedDomain( mesh, oStr ) << std::endl
         << "is (almost) zero" << std::endl;

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
  @param numberStoredFaces,
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
                       UInt &                 numberStoredFaces,
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
                   entityFlag_Type                  regionFlag,
                   bool                             verbose = false,
                   InternalEntitySelector           iSelect = InternalEntitySelector() )
{
    std::string line;

    Real x, y, z;

    UInt done = 0;
    UInt i;
    UInt numberVertices( 0 ), numberBoundaryVertices( 0 ),
         numberFaces   ( 0 ), numberBoundaryFaces   ( 0 ),
         numberPoints  ( 0 ), numberBoundaryPoints  ( 0 ),
         numberEdges   ( 0 ), numberBoundaryEdges   ( 0 );
    UInt numberVolumes ( 0 );
    UInt numberStoredFaces( 0 );
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
        std::cerr << " Error in readINRIAMeshFile = file " << fileName
                  << " not found or locked" << std::endl;
        std::abort();
    }

    if ( verbose )
    {
    std::cout << "Reading INRIA mesh file" << fileName << std::endl;
    }

    if ( ! readINRIAMeshFileHead( hstream, numberVertices, numberBoundaryVertices,
                                  numberBoundaryFaces, numberBoundaryEdges,
                                  numberVolumes, numberStoredFaces, shape, iSelect) )
    {
        std::cerr << " Error while reading INRIA mesh file headers" << std::endl;
        std::abort() ;
    }

    hstream.close();

    //Reopen the stream: I know it is stupid but this is how it goes
    std::ifstream myStream( fileName.c_str() );

    if ( myStream.fail() )
    {
        std::cerr << " Error in readINRIAMeshFile = file " << fileName
                  << " not found or locked" << std::endl;
        std::abort();
    }

    ASSERT_PRE0( GeoShape::S_shape == shape, "INRIA Mesh file and mesh element shape is not consistent" );

    // Euler formulas to get number of faces and number of edges
    numberFaces = 2 * numberVolumes + ( numberBoundaryFaces / 2 );
    Int num1  = numberVertices + numberVolumes;
    Int num2  = numberBoundaryVertices;
    Int num3  = numberBoundaryFaces;

    numberEdges = ( 3 * num3 - 2 * num2 ) / 4 + num1;

//    numberEdges = (int) numberVolumes + numberVertices + ( 3 * numberBoundaryFaces + dummy ) / 4;

    // Be a little verbose
    switch ( shape )
    {

    case HEXA:
        ASSERT_PRE0( GeoShape::S_numPoints == 8, "Sorry I can read only bilinear Hexa meshes" );
        std::cout << "Linear Hexa mesh" << std::endl;
        numberPoints =  numberVertices;
        numberBoundaryPoints = numberBoundaryVertices;
        break;

    case TETRA:
        if ( GeoShape::S_numPoints > 4 )
        {
            //if (GeoShape::S_numPoints ==6 )
            std::cout << "Quadratic Tetra mesh (from linear geometry)" << std::endl;
            numberPoints = numberVertices + numberEdges;

            /*
            numberBoundaryPoints=numberBoundaryVertices+numberBoundaryEdges;
            FALSE : numberBoundaryEdges is not known at this stage in a INRIA file
            (JFG 07/2002)
            I use the relation  numberBoundaryVertices + numberBoundaryFaces - 2 = numberBoundaryEdges,
            But, is it general (hole...)
            (JFG 07/2002)

            numberBoundaryEdges = ( Int ( numberBoundaryVertices + numberBoundaryFaces - Int (2  ) )
                                  > 0 ? ( numberBoundaryVertices + numberBoundaryFaces - 2 ) : 0 );
            numberBoundaryPoints = ( Int ( numberBoundaryVertices +
                                   ( numberBoundaryVertices + numberBoundaryFaces - Int ( 2 ) ) )
                                   > 0 ?numberBoundaryVertices +
                                   ( numberBoundaryVertices + numberBoundaryFaces - 2 ) : 0 );
             */
        }

        else
        {
            if ( verbose )
               {
                std::cout << "Linear Tetra Mesh" << std::endl;
               }

            numberPoints  = numberVertices;
            numberBoundaryPoints = numberBoundaryVertices;
            numberBoundaryEdges  = ( Int ( numberBoundaryVertices + numberBoundaryFaces - Int( 2 ) )
                                 > 0 ? ( numberBoundaryVertices + numberBoundaryFaces - 2 ) : 0 );
        }

        break;

    default:
        ERROR_MSG( "Current version of INRIA Mesh file reader only accepts TETRA and HEXA" );
    }

    oStr << "Number of Vertices        = "  << std::setw( 10 ) << numberVertices         << std::endl
         << "Number of BVertices       = "  << std::setw( 10 ) << numberBoundaryVertices << std::endl
         << "Number of Faces           = "  << std::setw( 10 ) << numberFaces            << std::endl
         << "Number of Boundary Faces  = "  << std::setw( 10 ) << numberBoundaryFaces    << std::endl
         << "Number of Stored Faces    = "  << std::setw( 10 ) << numberStoredFaces      << std::endl
         << "Number of Edges           = "  << std::setw( 10 ) << numberEdges            << std::endl
         << "Number of Boundary Edges  = "  << std::setw( 10 ) << numberBoundaryEdges    << std::endl
         << "Number of Points          = "  << std::setw( 10 ) << numberPoints           << std::endl
         << "Number of Boundary Points = "  << std::setw( 10 ) << numberBoundaryPoints   << std::endl
         << "Number of Volumes         = "  << std::setw( 10 ) << numberVolumes          << std::endl;

    // Set all basic data structure

    // I store all Points
    mesh.setMaxNumPoints       ( numberPoints, true );
    mesh.setMaxNumGlobalPoints ( numberPoints );
    mesh.setNumBPoints         ( numberBoundaryPoints );
    mesh.setNumVertices        ( numberVertices );
    mesh.setNumGlobalVertices  ( numberVertices );
    mesh.setNumBVertices       ( numberBoundaryVertices );
    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges        ( numberBoundaryEdges );
    mesh.setMaxNumGlobalEdges  ( numberEdges );
    mesh.setNumEdges           ( numberEdges ); // Here the REAL number of edges (all of them)
    mesh.setNumBEdges          ( numberBoundaryEdges );
    // Only Boundary Faces
    mesh.setMaxNumFaces        ( numberStoredFaces );
    mesh.setMaxNumGlobalFaces  ( numberBoundaryFaces );
    mesh.setNumFaces           ( numberFaces ); // Here the REAL number of faces (all of them)
    mesh.setNumBFaces          ( numberBoundaryFaces );

    mesh.setMaxNumVolumes      ( numberVolumes, true );
    mesh.setMaxNumGlobalVolumes( numberVolumes );

    mesh.setMarker             ( regionFlag ); // Add Marker to list of Markers

    typedef typename RegionMesh3D<GeoShape, MC>::PointType  PointType;
    typedef typename RegionMesh3D<GeoShape, MC>::VolumeType VolumeType;


    typename RegionMesh3D<GeoShape, MC>::PointType  * pointerPoint  = 0;
    typename RegionMesh3D<GeoShape, MC>::EdgeType   * pointerEdge   = 0;
    typename RegionMesh3D<GeoShape, MC>::FaceType   * pointerFace   = 0;
    typename RegionMesh3D<GeoShape, MC>::VolumeType * pointerVolume = 0;
    // addPoint()/Face()/Edge() returns a reference to the last stored point
    // I use that information to set all point info, by using a pointer.

    UInt count = 0;
    Int  ibc;

    // To account for internal faces
    if ( numberStoredFaces > numberBoundaryFaces )
    {
        faceHelp.resize( numberStoredFaces - numberBoundaryFaces );
        faceHelpIterator = faceHelp.begin();

        oStr << "WARNING: The mesh file (apparently) contains "
             << numberStoredFaces - numberBoundaryFaces << " internal faces" << std::endl;

    }

    while ( nextGoodLine( myStream, line ).good() )
    {
        if ( line.find( "Vertices" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), myStream );

            for ( i = 0; i < numberVertices; i++ )
            {
                myStream >> x >> y >> z >> ibc;

//                if (ibc == 1 ) ibc = 100;

                if ( !iSelect(entityFlag_Type(ibc)))
                {
                    ++count;
                // Boundary point. Boundary switch set by the mesh method.
                    pointerPoint = &mesh.addPoint( true );
                    pointerPoint->setMarker( entityFlag_Type( ibc ) );
                }

                else
                {
                    pointerPoint = &mesh.addPoint( false );
                }

                pointerPoint->setId     ( i + 1 );
                pointerPoint->setLocalId( i + 1 );
                pointerPoint->x() = x;
                pointerPoint->y() = y;
                pointerPoint->z() = z;
                pointerPoint->setMarker( entityFlag_Type( ibc ) );

                mesh.localToGlobalNode().insert( std::make_pair( i + 1, i + 1 ) );
                mesh.globalToLocalNode().insert( std::make_pair( i + 1, i + 1 ) );
            }

            oStr << "Vertices read " << std::endl;
            oStr << "Size of the node storage is "
                 << count * sizeof( PointType ) / 1024. / 1024. << std::endl;
            done++;

            if ( count != numberBoundaryVertices )
              {
                std::cerr << "Number boundary points inconsistent!" << std::endl;
              }
        }

        if ( line.find( "Triangles" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), myStream );
            oStr << "Reading boundary faces " << std::endl;

            for ( i = 0; i < numberStoredFaces; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> ibc;

                if ( numberStoredFaces > numberBoundaryFaces )
                {
                    if ( mesh.point( p1 ).boundary() && mesh.point( p2 ).boundary() &&
                            mesh.point( p3 ).boundary() )
                    {
                        pointerFace = &( mesh.addFace( true ) ); // Boundary faces
                        pointerFace->setMarker( entityFlag_Type( ibc ) );
                        pointerFace->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                        pointerFace->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                        pointerFace->setPoint( 3, mesh.point( p3 ) ); // set face conn.

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

                    pointerFace = &( mesh.addFace( true ) ); // Only boundary faces

                    pointerFace->setMarker( entityFlag_Type( ibc ) );
                    pointerFace->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                    pointerFace->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                    pointerFace->setPoint( 3, mesh.point( p3 ) ); // set face conn.
                }
            }

         for ( faceHelpIterator = faceHelp.begin();
               faceHelpIterator != faceHelp.end(); ++faceHelpIterator )
          {
              p1  = faceHelpIterator->i1;
              p2  = faceHelpIterator->i2;
              p3  = faceHelpIterator->i3;
              ibc = faceHelpIterator->ibc;
              pointerFace  = &( mesh.addFace( false ) ); // INTERNAL FACE
              pointerFace->setMarker( entityFlag_Type( ibc ) );
              pointerFace->setPoint( 1, mesh.point( p1 ) ); // set face conn.
              pointerFace->setPoint( 2, mesh.point( p2 ) ); // set face conn.
              pointerFace->setPoint( 3, mesh.point( p3 ) ); // set face conn.
          }

            oStr << "Boundary faces read " << std::endl;
            done++;
        }

        if ( line.find( "Quadrilaterals" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), myStream );

            oStr << "Reading boundary faces " << std::endl;

            for (UInt i = 0; i < numberBoundaryFaces; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> p4 >> ibc;

                if ( numberStoredFaces > numberBoundaryFaces )
                {
                    if ( mesh.point( p1 ).boundary() && mesh.point( p2 ).boundary() &&
                            mesh.point( p3 ).boundary() )
                    {
                        pointerFace = &( mesh.addFace( true ) ); // Boundary faces
                        pointerFace->setMarker( entityFlag_Type( ibc ) );
                        pointerFace->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                        pointerFace->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                        pointerFace->setPoint( 3, mesh.point( p3 ) ); // set face conn.
                        pointerFace->setPoint( 4, mesh.point( p4 ) ); // set face conn.

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
                    pointerFace = &( mesh.addFace( true ) ); // Only boundary faces
                    pointerFace->setMarker( entityFlag_Type( ibc ) );
                    pointerFace->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                    pointerFace->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                    pointerFace->setPoint( 3, mesh.point( p3 ) ); // set face conn.
                    pointerFace->setPoint( 4, mesh.point( p4 ) ); // set face conn.
                }
            }

            oStr << "Boundary faces read " << std::endl;

            for ( faceHelpIterator=faceHelp.begin();
                  faceHelpIterator!=faceHelp.end(); ++faceHelpIterator )
            {
                p1 = faceHelpIterator->i1;
                p2 = faceHelpIterator->i2;
                p3 = faceHelpIterator->i3;
                p4 = faceHelpIterator->i4;
                ibc = faceHelpIterator->ibc;
                pointerFace = &( mesh.addFace( false ) ); // INTERNAL FACE
                pointerFace->setMarker( entityFlag_Type( ibc ) );
                pointerFace->setPoint( 1, mesh.point( p1 ) ); // set face conn.
                pointerFace->setPoint( 2, mesh.point( p2 ) ); // set face conn.
                pointerFace->setPoint( 3, mesh.point( p3 ) ); // set face conn.
                pointerFace->setPoint( 4, mesh.point( p4 ) ); // set face conn.
            }
            done++;
        }

        if ( line.find( "Edges" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), myStream );
            oStr << "Reading boundary edges " << std::endl;

            for ( i = 0; i < numberBoundaryEdges; i++ )
            {
                myStream >> p1 >> p2 >> ibc;
                pointerEdge = &mesh.addEdge( true ); // Only boundary edges.
                pointerEdge->setMarker( entityFlag_Type( ibc ) );
                pointerEdge->setPoint( 1, mesh.point( p1 ) ); // set edge conn.
                pointerEdge->setPoint( 2, mesh.point( p2 ) ); // set edge conn.
            }
            oStr << "Boundary edges read " << std::endl;
            done++;
        }

        if ( line.find( "Tetrahedra" ) != std::string::npos )
        {
            count = 0;
            nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), myStream );
            oStr << "Reading volumes " << std::endl;

            for ( i = 0; i < numberVolumes; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> p4 >> ibc;
                pointerVolume = &mesh.addVolume();
                pointerVolume->setId     ( i + 1 );
                pointerVolume->setLocalId( i + 1);
                pointerVolume->setPoint( 1, mesh.point( p1 ) );
                pointerVolume->setPoint( 2, mesh.point( p2 ) );
                pointerVolume->setPoint( 3, mesh.point( p3 ) );
                pointerVolume->setPoint( 4, mesh.point( p4 ) );
                pointerVolume->setMarker( entityFlag_Type( ibc ) );
//                mesh.localToGlobalElem().insert(std::make_pair(i+1, i+1));
//                mesh.globalToLocalElem().insert(std::make_pair(i+1, i+1));
                count++;
            }
            oStr << "size of the volume storage is " << sizeof( VolumeType ) * count / 1024. / 1024.
                 << " Mo." << std::endl;
            oStr << count << " Volume elements read" << std::endl;
            done++;
        }

        if ( line.find( "Hexahedra" ) != std::string::npos )
        {
            count = 0;
            nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), myStream );
            oStr << "Reading volumes " << std::endl;
            for ( i = 0; i < numberVolumes; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> p8 >> ibc;
                pointerVolume = &mesh.addVolume();
                pointerVolume->setId     ( i + 1 );
                pointerVolume->setLocalId( i + 1);
//                pointerVolume->id() = i + 1;
                pointerVolume->setPoint( 1, mesh.point( p1 ) );
                pointerVolume->setPoint( 2, mesh.point( p2 ) );
                pointerVolume->setPoint( 3, mesh.point( p3 ) );
                pointerVolume->setPoint( 4, mesh.point( p4 ) );
                pointerVolume->setPoint( 5, mesh.point( p5 ) );
                pointerVolume->setPoint( 6, mesh.point( p6 ) );
                pointerVolume->setPoint( 7, mesh.point( p7 ) );
                pointerVolume->setPoint( 8, mesh.point( p8 ) );
                pointerVolume->setMarker( entityFlag_Type( ibc ) );

                count++;
            }
            oStr << count << " Volume elements read" << std::endl;
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
    if ( shape == TETRA && GeoShape::S_numPoints > 4 )
    {
        MeshUtility::p2MeshFromP1Data( mesh );
    }

    myStream.close();

    Real vols[ 3 ];
    getVolumeFromFaces( mesh, vols, oStr );

    oStr << "Volume enclosed by the mesh computed by integration on boundary faces"  << std::endl;
    oStr << "INT(X)     INT(Y)      INT(Z)     <- they should be equal and equal to" << std::endl
         << "                                     the voulume enclosed by the mesh"  << std::endl;
    oStr << vols[ 0 ] << "     " << vols[ 1 ] << "     " << vols[ 2 ] << std::endl;
    oStr << "Boundary faces are defining a closed surface if "
         << testClosedDomain( mesh, oStr ) << std::endl
         << "is (almost) zero" << std::endl;

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
              entityFlag_Type              regionFlag,
              bool                         verbose = false )
{
    std::ifstream inputFile ( fileName.c_str() );

#ifdef DEBUG
    std::debug( 3000 ) << "Gmsh reading: " << fileName << "\n";
#endif

    //    char buffer[256];
    std::string buffer;

    for (Int ii = 0; ii < 6; ++ii)
    {
        inputFile >> buffer;
        std::cout << "buffer = " << buffer << "\n";
    }

    UInt numberNodes;
    inputFile >> numberNodes;

#ifdef DEBUG
    std::debug( 3000 ) << "Number of nodes = " << numberNodes;
#endif

    // Add Marker to list of Markers
    mesh.setMarker( regionFlag );


    std::vector<Real> x( 3 * numberNodes );
    std::vector<bool> isonboundary( numberNodes );
    std::vector<UInt> whichboundary( numberNodes );

#ifdef DEBUG
    std::debug( 3000 ) << "Reading " << numberNodes << " nodes\n";
#endif

    std::map<Int,Int> itoii;

    for ( UInt i = 0; i < numberNodes; ++i )
    {
        UInt ni;
        inputFile >> ni
        >> x[ 3 * i ]
        >> x[ 3 * i + 1 ]
        >> x[ 3 * i + 2 ];

        itoii[ ni - 1] = i;
    }
    inputFile >> buffer;

#ifdef DEBUG
    std::debug( 3000 ) << "buffer = " << buffer << "\n";
#endif

    inputFile >> buffer;

#ifdef DEBUG
    std::debug( 3000 ) << "buffer = " << buffer << "\n";
#endif

    UInt numberElements;
    inputFile >> numberElements;

    typename RegionMesh3D<GeoShape, MC>::EdgeType   * pointerEdge   = 0;
    typename RegionMesh3D<GeoShape, MC>::FaceType   * pointerFace   = 0;
    typename RegionMesh3D<GeoShape, MC>::VolumeType * pointerVolume = 0;

#ifdef DEBUG
    std::debug( 3000 ) << "number of elements: " << numberElements << "\n";
#endif

    std::vector<std::vector<int> > e( numberElements );
    std::vector<int>               et( numberElements );
    std::vector<int>               etype( numberElements );
    std::vector<int>               gt( 32 );
    gt.assign( 32, 0 );

    for ( UInt i = 0; i < numberElements; ++i )
    {
        Int ne, t, np;

        //std::debug() << i + 1 << " ";

        inputFile >> buffer;
        inputFile >> ne;

        //std::debug() << ne << " ";

        switch ( ne )
        {
        case(2):
            np = 3;
            break;

        case(4):
            np = 4;
            break;

        case(15):
            np = 1;
            break;

        default:
            np = 0;

#ifdef DEBUG
            std::debug( 3000 ) << "Element type unknown " << ne << "\n";
#endif

            ASSERT( true, "Elements type unsupported.\n" )
        }

        inputFile >> t;

        //std::debug() << t << " ";

        bool ibcSet = false;
        Int  flag   = 0;
        Int tag( 0 );

        for ( Int iflag = 0; iflag < t; ++iflag )
        {
            inputFile >> flag;

            if ( !ibcSet )
            {
                tag = flag;
                ibcSet = true;
            }
        }

        ++gt[ ne ];

        etype[ i ] = ne;

        et[ i ] = tag;
        e[ i ].resize( np );

        Int p = 0;

        while ( p != np )
        {
            Int node;
            inputFile >> node;
            e[ i ][ p ] = node;
            e[ i ][ p ] = itoii[ e[ i ][ p ] - 1];
            e[ i ][ p ] += 1;

            ++p;
        }
    }


    // Euler formulas
    UInt n_volumes = gt[ 4 ];
    UInt n_faces_boundary = gt[ 2 ];
    UInt n_faces_total = 2 * n_volumes + ( n_faces_boundary / 2 );

    mesh.setMaxNumGlobalPoints( numberNodes );
    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumEdges( gt[ 1 ] );
    mesh.setNumEdges   ( gt[ 1 ] ); // Here the REAL number of edges (all of them)
    mesh.setNumBEdges  ( gt[ 1 ] );
    mesh.setMaxNumGlobalEdges( gt[ 1 ] );

#ifdef DEBUG
    std::debug( 3000 ) << "number of edges= " << gt[ 1 ] << "\n";
#endif

    // Only Boundary Faces
    mesh.setMaxNumFaces( n_faces_total );
    mesh.setNumFaces   ( n_faces_total ); // Here the REAL number of edges (all of them)
    //mesh.setMaxNumFaces( n_faces_boundary );
    //mesh.setNumFaces   ( n_faces_boundary ); // Here the REAL number of edges (all of them)
    mesh.setNumBFaces  ( n_faces_boundary );
    mesh.setMaxNumGlobalFaces( n_faces_total );

#ifdef DEBUG
    std::debug( 3000 ) << "number of faces = " << n_faces_boundary << "\n";
#endif

    mesh.setMaxNumVolumes( n_volumes, true );
    mesh.setMaxNumGlobalVolumes( n_volumes );

#ifdef DEBUG
    std::debug( 3000 ) << "number of volumes = " << n_volumes << "\n";
#endif

    isonboundary.assign( numberNodes, false );
    whichboundary.assign( numberNodes, 0 );

    for ( UInt i = 0; i < numberElements; ++i )
    {
        switch ( etype[ i ] )
        {
            // triangular faces (linear)
        case 2:
        {
            isonboundary[ e[ i ][ 0 ] - 1 ] = true;
            isonboundary[ e[ i ][ 1 ] - 1 ] = true;
            isonboundary[ e[ i ][ 2 ] - 1 ] = true;

            whichboundary[ e[ i ][ 0 ] - 1 ] = et[ i ];
            whichboundary[ e[ i ][ 1 ] - 1 ] = et[ i ];
            whichboundary[ e[ i ][ 2 ] - 1 ] = et[ i ];
        }
        }
    }
    // add the point to the mesh
    typename RegionMesh3D<GeoShape, MC>::PointType * pointerPoint = 0;

    mesh.setMaxNumPoints( numberNodes, true );
    mesh.setNumVertices ( numberNodes );
    mesh.setNumBVertices( std::count( isonboundary.begin(), isonboundary.end(), true ) );
    mesh.setNumBPoints  ( mesh.numBVertices() );

#ifdef DEBUG
    std::debug( 3000 ) << "number of points : "            << mesh.numPoints() << "\n";
    std::debug( 3000 ) << "number of boundary points : "   << mesh.numBPoints() << "\n";
    std::debug( 3000 ) << "number of vertices : "          << mesh.numVertices() << "\n";
    std::debug( 3000 ) << "number of boundary vertices : " << mesh.numBVertices() << "\n";
#endif

    for ( UInt i = 0; i < numberNodes; ++i )
    {
        pointerPoint = &mesh.addPoint( isonboundary[ i ] );
        pointerPoint->setMarker( whichboundary[ i ] );
        pointerPoint->setId     ( i + 1 );
        pointerPoint->setLocalId( i + 1 );
        pointerPoint->x() = x[ 3 * i ];
        pointerPoint->y() = x[ 3 * i + 1 ];
        pointerPoint->z() = x[ 3 * i + 2 ];
        mesh.localToGlobalNode().insert( std::make_pair( i + 1, i + 1 ) );
        mesh.globalToLocalNode().insert( std::make_pair( i + 1, i + 1 ) );
    }

    Int numberVolumes = 1;
    // add the element to the mesh
    for ( UInt i = 0; i < numberElements; ++i )
    {
        switch ( etype[ i ] )
        {
        // segment(linear)
        case 1:
        {
            pointerEdge = &( mesh.addEdge( true ) );
            pointerEdge->setMarker( entityFlag_Type( et[ i ] ) );
            pointerEdge->setPoint( 1, mesh.point( e[ i ][ 0 ] ) );
            pointerEdge->setPoint( 2, mesh.point( e[ i ][ 1 ] ) );



        }
        break;

        // triangular faces (linear)
        case 2:
        {
            pointerFace = &( mesh.addFace( true ) );
            pointerFace->setMarker( entityFlag_Type( et[ i ] ) );
            pointerFace->setPoint( 1, mesh.point( e[ i ][ 0 ] ) );
            pointerFace->setPoint( 2, mesh.point( e[ i ][ 1 ] ) );
            pointerFace->setPoint( 3, mesh.point( e[ i ][ 2 ] ) );

        }
        break;

        // quadrangular faces(linear)
        case 3:
        {
            pointerFace = &( mesh.addFace( true ) );
            pointerFace->setMarker( entityFlag_Type( et[ i ] ) );
            pointerFace->setPoint( 1, mesh.point( e[ i ][ 0 ] ) );
            pointerFace->setPoint( 2, mesh.point( e[ i ][ 1 ] ) );
            pointerFace->setPoint( 3, mesh.point( e[ i ][ 2 ] ) );
            pointerFace->setPoint( 4, mesh.point( e[ i ][ 3 ] ) );
        }
        break;

        // tetrahedrons(linear)
        case 4:
        {
            pointerVolume = &( mesh.addVolume() );
            pointerVolume->setId     ( numberVolumes );
            pointerVolume->setLocalId( numberVolumes++ );
//                pointerVolume->id() = i + 1;
            pointerVolume->setMarker( entityFlag_Type( et[ i ] ) );
            pointerVolume->setPoint( 1, mesh.point( e[ i ][ 0 ] ) );
            pointerVolume->setPoint( 2, mesh.point( e[ i ][ 1 ] ) );
            pointerVolume->setPoint( 3, mesh.point( e[ i ][ 2 ] ) );
            pointerVolume->setPoint( 4, mesh.point( e[ i ][ 3 ] ) );
        }
        break;

        // hexahedrons(linear)
        case 5:
        {
            pointerVolume = &( mesh.addVolume() );

            pointerVolume->setId     ( i + 1 );
            pointerVolume->setLocalId( i + 1 );
//                pointerVolume->id() = i + 1;
            pointerVolume->setMarker( entityFlag_Type( et[ i ] ) );
            pointerVolume->setPoint( 1, mesh.point( e[ i ][ 0 ] ) );
            pointerVolume->setPoint( 2, mesh.point( e[ i ][ 1 ] ) );
            pointerVolume->setPoint( 3, mesh.point( e[ i ][ 2 ] ) );
            pointerVolume->setPoint( 4, mesh.point( e[ i ][ 3 ] ) );
            pointerVolume->setPoint( 5, mesh.point( e[ i ][ 4 ] ) );
            pointerVolume->setPoint( 6, mesh.point( e[ i ][ 5 ] ) );
            pointerVolume->setPoint( 7, mesh.point( e[ i ][ 6 ] ) );
            pointerVolume->setPoint( 8, mesh.point( e[ i ][ 7 ] ) );
        }
        break;
        }
    }

    Switch sw;

    std::stringstream discardedLog;
    std::ostream& oStr = verbose ? std::cout : discardedLog;

    if ( checkMesh3D(mesh, sw, true,verbose,oStr,std::cerr,oStr) == false )
    {
        std::ostringstream ex;
        ex << "invalid mesh from GSMH";

        throw std::logic_error( ex.str() );
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
               entityFlag_Type             regionFlag,
               bool                        verbose = false )
{
    // I will extract lines from iostream
    std::string line;

    // number of Geo Elements
    UInt numberVertices( 0 ), numberBoundaryVertices( 0 ),
         numberPoints  ( 0 ), numberBoundaryPoints  ( 0 ),
         numberEdges   ( 0 ), numberBoundaryEdges   ( 0 ),
         numberFaces   ( 0 ), numberBoundaryFaces   ( 0 ),
         numberVolumes ( 0 );

    // During the first access to file, build list of structures
    Vector pointCoordinates;
    std::vector<UInt> edgePointID, facePointID, volumePointID;

    // for poit i-th: is it a boundary point?
    std::vector<bool> boundaryPoint;

    // build a list of boundary edges, since netgen is not writing all of them
    BareItemsHandler<BareEdge> bihBedges;

    // flags for boundary entities
    std::vector<entityFlag_Type> bcnsurf, bcnpoints;

    // bitstream to check which file section has already been visited
    UInt flag;

    typename MC::PointMarker pointMarker;
    typename MC::EdgeMarker edgeMarker;

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
            std::stringstream parseLine( line );

            parseLine >> numberVertices;
            std::cout << "[readNetgenMesh] found " << numberVertices << " vertices... " << std::flush;

            bcnpoints.resize( numberVertices + 1 );
            boundaryPoint.resize( numberVertices + 1 );
            pointCoordinates.resize( nDimensions * numberVertices );

            for ( UInt i = 0; i < numberVertices; i++ )
            {
                // helping variables to read netgen file fields
                Real x, y, z;

                getline( fstreamp, line );
                std::stringstream parseLine( line );

                parseLine >> x >> y >> z;

                pointCoordinates[ i * nDimensions ]     = x;
                pointCoordinates[ i * nDimensions + 1 ] = y;
                pointCoordinates[ i * nDimensions + 2 ] = z;
                boundaryPoint  [ i ] = false;
                bcnpoints[ i ] = S_NULLFLAG;
            }
            boundaryPoint  [ numberVertices ] = false;
            bcnpoints[ numberVertices ] = S_NULLFLAG;

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
            std::stringstream parseLine( line );

            parseLine >> numberBoundaryEdges; // this will be discarded in 3D meshes

            std::cout << "[readNetgenMesh] found " << numberBoundaryEdges << " boundary edges... " << std::flush;

            edgePointID.resize( 2 * numberBoundaryEdges );

            for ( UInt i = 0; i < numberBoundaryEdges; ++i)
            {
                UInt surfnr, t, p1, p2;

                getline( fstreame, line );
                std::stringstream parseLine( line );

                parseLine >> surfnr >> t >> p1 >> p2;

                edgePointID[ 2 * i ] = p1;
                edgePointID[ 2 * i + 1 ] = p2;

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
            std::stringstream parseLine( line );

            parseLine >> numberVolumes;
            std::cout << "[readNetgenMesh] found " << numberVolumes << " volumes... " << std::flush;

            volumePointID.resize( 4 * numberVolumes );

            for ( UInt i = 0; i < numberVolumes; i++ )
            {
                UInt matnr, np, p1, p2, p3, p4;

                getline( fstreamv, line );
                std::stringstream parseLine( line );

                parseLine >> matnr >> np >> p1 >> p2 >> p3 >> p4;

                volumePointID[ 4 * i ] = p1;
                volumePointID[ 4 * i + 1 ] = p2;
                volumePointID[ 4 * i + 2 ] = p3;
                volumePointID[ 4 * i + 3 ] = p4;

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
            std::stringstream parseLine( line );

            parseLine >> numberBoundaryFaces;

            std::cout << "[readNetgenMesh] found " << numberBoundaryFaces
                      << " boundary faces... "     << std::flush;

            facePointID.resize( 3 * numberBoundaryFaces );
            bcnsurf.resize( numberBoundaryFaces + 1 );
            bcnsurf[ 0 ] = S_NULLFLAG;


            for ( UInt i = 0; i < numberBoundaryFaces; i++ )
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
                std::stringstream parseLine( line );

                parseLine >> surfnr >> bcnr >> domin >> domout >> np >> p1 >> p2 >> p3;

                facePointID[ 3 * i ] = p1;
                facePointID[ 3 * i + 1 ] = p2;
                facePointID[ 3 * i + 2 ] = p3;

                bcnsurf[ i + 1 ] = bcnr;
                //          std::cout<<"[readNetgenMesh] bcnr = " << bcnr << std::endl;
                ASSERT( np==3, "Error in readNetgenMesh: only triangular surfaces supported" )

                //assume p1!=p2!=p3!=p1
                numberBoundaryVertices += boundaryPoint[ p1 ] ? 0:1;
                numberBoundaryVertices += boundaryPoint[ p2 ] ? 0:1;
                numberBoundaryVertices += boundaryPoint[ p3 ] ? 0:1;
                boundaryPoint[ p1 ] = boundaryPoint[ p2 ] = boundaryPoint[ p3 ] = true;

                /* here I set the boundary points marker
                   note: this works only with my patch
                   of strongerFlag
                   Face flag is assigned to face points
                   A point receives the "stronger" flag of the faces it belongs to
                    */
                bcnpoints[ p1 ] = pointMarker.setStrongerMarker( bcnpoints[ p1 ], bcnr );
                bcnpoints[ p2 ] = pointMarker.setStrongerMarker( bcnpoints[ p2 ], bcnr );
                bcnpoints[ p3 ] = pointMarker.setStrongerMarker( bcnpoints[ p3 ], bcnr );

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

                bihBedges.addIfNotThere( bed, ( ID )S_NULLFLAG );
                bihBedges[ bed ] = ( ID )edgeMarker.setStrongerMarker( bihBedges[ bed ], bcnr );

                bed = setBareEdge( p2, p3 );
                bihBedges.addIfNotThere( bed,( ID )S_NULLFLAG );
                bihBedges[ bed ] = ( ID )edgeMarker.setStrongerMarker( bihBedges[ bed ], bcnr );

                bed = setBareEdge( p3, p1 );
                bihBedges.addIfNotThere( bed,( ID )S_NULLFLAG );
                bihBedges[ bed ] = ( ID )edgeMarker.setStrongerMarker( bihBedges[ bed ], bcnr );
            }
            flag&=~4;
            // in the 3D case the only way to know the number of edges on boundary faces
            // is to count them!
            numberBoundaryEdges = bihBedges.howMany();
            std::cout << "loaded." << std::endl;
            break;
        }
    }

    fstreamf.close();

    ASSERT( flag==0, "[readNetgenMesh] the mesh file does not have all the required sections." )

    std::cout << "[readNetgenMesh] computed " << numberBoundaryVertices << " boundary vertices" << std::endl;

    // Euler formulas
    numberFaces = 2 * numberVolumes + ( numberBoundaryFaces / 2 );
    numberEdges = numberVolumes + numberVertices + ( 3 * numberBoundaryFaces - 2 * numberBoundaryVertices ) / 4;

    // Be a little verbose
    if ( GeoShape::S_numPoints > 4 )
    {
        std::cout << "Quadratic Tetra  Mesh (from Linear geometry)" <<std::endl;
        numberPoints = numberVertices + numberEdges;
        numberBoundaryPoints = numberBoundaryVertices + numberBoundaryEdges;    // I calculated the real numberBoundaryEdges before
    }

    else
    {
        std::cout << "Linear Tetra Mesh" <<std::endl;
        numberPoints = numberVertices;
        numberBoundaryPoints = numberBoundaryVertices;
    }

    std::stringstream discardedLog;
    std::ostream& oStr = verbose ? std::cout : discardedLog;
    //points can be only vertices or on edges too

    std::cout << "Number of Vertices        = " << std::setw( 10 ) << numberVertices  << std::endl
              << "Number of BVertices       = " << std::setw( 10 ) << numberBoundaryVertices << std::endl;
    oStr      << "Number of Faces           = " << std::setw( 10 ) << numberFaces  << std::endl
              << "Number of Boundary Faces  = " << std::setw( 10 ) << numberBoundaryFaces << std::endl
              << "Number of Edges           = " << std::setw( 10 ) << numberEdges  << std::endl
              << "Number of Boundary Edges  = " << std::setw( 10 ) << numberBoundaryEdges << std::endl;
    std::cout << "Number of Points          = " << std::setw( 10 ) << numberPoints  << std::endl
              << "Number of Boundary Points = " << std::setw( 10 ) << numberBoundaryPoints << std::endl
              << "Number of Volumes         = " << std::setw( 10 ) << numberVolumes  << std::endl;

    // Set all basic data structure

    // I store all Points
    mesh.setMaxNumPoints       ( numberPoints, true );
    mesh.setMaxNumGlobalPoints ( numberPoints );
    mesh.setNumBPoints         ( numberBoundaryPoints );
    mesh.setNumVertices        ( numberVertices );
    mesh.setNumGlobalVertices  ( numberVertices );
    mesh.setNumBVertices       ( numberBoundaryVertices );

    // Only Boundary Edges
    mesh.setMaxNumEdges        ( numberBoundaryEdges );
    mesh.setMaxNumGlobalEdges  ( numberBoundaryEdges );
    mesh.setNumEdges           ( numberEdges ); // Here the REAL number of edges (all of them)
    mesh.setNumBEdges          ( numberBoundaryEdges );

    // Only Boundary Faces
    mesh.setMaxNumFaces        ( numberBoundaryFaces );
    mesh.setMaxNumGlobalFaces  ( numberBoundaryFaces );
    mesh.setNumFaces           ( numberFaces ); // Here the REAL number of faces (all of them)
    mesh.setNumBFaces          ( numberBoundaryFaces );

    mesh.setMaxNumVolumes      ( numberVolumes, true );
    mesh.setMaxNumGlobalVolumes( numberVolumes );

    mesh.setMarker             ( regionFlag ); // Add Marker to list of Markers

    typename RegionMesh3D<GeoShape,MC>::PointType  * pointerPoint  = 0;
    typename RegionMesh3D<GeoShape,MC>::EdgeType   * pointerEdge   = 0;
    typename RegionMesh3D<GeoShape,MC>::FaceType   * pointerFace   = 0;
    typename RegionMesh3D<GeoShape,MC>::VolumeType * pointerVolume = 0;

    // addPoint()/Face()/Edge() returns a reference to the last stored point
    // I use that information to set all point info, by using a pointer.

    std::cout << "[readmesh3D] boundaryPoint.size() = " << boundaryPoint.size()
              << ", bcnpoints.size() = "          << bcnpoints.size()
              << std::endl;

    for ( UInt i = 0; i < numberVertices; i++ )
    {
        pointerPoint=&mesh.addPoint( boundaryPoint[ i + 1 ] ); //true if boundary point

        pointerPoint->setId     ( i + 1 );
        pointerPoint->setLocalId( i + 1 );
        mesh.localToGlobalNode().insert( std::make_pair( i + 1, i + 1 ) );
        mesh.globalToLocalNode().insert( std::make_pair( i + 1, i + 1 ) );

        pointerPoint->setMarker( bcnpoints[ i + 1 ] );
        pointerPoint->x() = pointCoordinates[ nDimensions * i ];
        pointerPoint->y() = pointCoordinates[ nDimensions * i + 1 ];
        pointerPoint->z() = pointCoordinates[ nDimensions * i + 2 ];
    }

    std::cout << "[readmesh3D] added points." << std::endl;

    /*
       here I set the real boundary edges that I stored
       in bihBedges
    */

    BareItemsHandler<BareEdge>::const_iterator bedge = bihBedges.begin();

    for ( UInt i=0; i < numberBoundaryEdges; i++ )
    {
        UInt p1, p2;

        pointerEdge = &mesh.addEdge( true ); // Only boundary edges.
        pointerEdge->setMarker( entityFlag_Type( bedge->second ) );
        p1 = bedge->first.first;
        p2 = bedge->first.second;
        pointerEdge->setPoint( 1, mesh.point( p1 ) ); // set edge conn.
        pointerEdge->setPoint( 2, mesh.point( p2 ) ); // set edge conn.

        bedge++;
    }

    std::cout << "[readmesh3D] added edges." << std::endl;

    for ( UInt i = 0; i < numberVolumes; i++ )
    {
        UInt p1, p2, p3, p4;

        pointerVolume = &mesh.addVolume();
        pointerVolume->setId( i + 1 );
        pointerVolume->setLocalId(i + 1 );
        p1 = volumePointID[ 4 * i ];
        p2 = volumePointID[ 4 * i + 1 ];
        p3 = volumePointID[ 4 * i + 2 ];
        p4 = volumePointID[ 4 * i + 3 ];
        pointerVolume->setPoint( 1, mesh.point( p1 ) );
        pointerVolume->setPoint( 2, mesh.point( p2 ) );
        pointerVolume->setPoint( 3, mesh.point( p3 ) );
        pointerVolume->setPoint( 4, mesh.point( p4 ) );
    }

    std::cout << "[readmesh3D] added volumes." << std::endl;

    for ( UInt i = 0; i < numberBoundaryFaces; i++ )
    {
        UInt p1, p2, p3;

        pointerFace = &mesh.addFace( true ); // Only boundary faces
        p1 = facePointID[ 3 * i ];
        p2 = facePointID[ 3 * i + 1 ];
        p3 = facePointID[ 3 * i + 2 ];

        pointerFace->setMarker( entityFlag_Type( bcnsurf[ i + 1 ] ) );
        pointerFace->setPoint( 1, mesh.point( p1 ) ); // set face conn.
        pointerFace->setPoint( 2, mesh.point( p2 ) ); // set face conn.
        pointerFace->setPoint( 3, mesh.point( p3 ) ); // set face conn.
    }
    std::cout << "[readmesh3D] added faces." << std::endl;

    // This part is to build a P2 mesh from a P1 geometry

    if ( GeoShape::S_numPoints > 4 )
    {
    MeshUtility::p2MeshFromP1Data( mesh );
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

    getVolumeFromFaces( mesh, vols, oStr );

    oStr << "Volume enclosed by the mesh computed by integration on boundary faces" << std::endl
         << "INT(X)     INT(Y)      INT(Z)      <- they should be equal and equal to" << std::endl
         << "                                   the volume enclosed by the mesh " << std::endl
         << vols[ 0 ] << "      " << vols[ 1 ] << "      " << vols[ 2 ] << std::endl
         << "Boundary faces are defining aclosed surface if "
         << testClosedDomain( mesh, oStr ) << std::endl
         << " is (almost) zero" << std::endl;

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
