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
    @brief new implementation of 3D mesh readers

    @author Antonio Cervone <ant.cervone@gmail.com>

    @date 09-08-2012
*/

#ifndef _IMPORTERMESH3DNEW_HH_
#define _IMPORTERMESH3DNEW_HH_ 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/MeshElementBare.hpp>
#include <lifev/core/mesh/MeshChecks.hpp>
#include <lifev/core/mesh/InternalEntitySelector.hpp>
#include <lifev/core/mesh/BareMesh.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/filter/ImporterMesh3D.hpp>

namespace LifeV
{

//! convertBareMesh - convert a previously read BareMesh in a RegionMesh object
/*!
  Starting from a BareMesh, this routine generates a fully compliant RegionMesh object

  @param bareMesh, the bare mesh data structure in input.
  @param mesh, the mesh data structure to fill in.
  @param regionFlag, the identifier for the region.
  @param verbose, setting it as true, the output is verbose (the default is false).
  @param iSelect,
  @return true if everything went fine, false otherwise.
*/

template <typename GeoShape, typename MC>
bool
convertBareMesh ( RegionMeshBare<GeoShape> &  bareMesh,
                  RegionMesh<GeoShape, MC>& mesh,
                  bool                        verbose = false,
                  InternalEntitySelector      iSelect = InternalEntitySelector() )
{
    Real x, y, z;
    UInt p[ GeoShape::S_numVertices ];

    UInt done = 0;
    UInt i;
    UInt numberVertices         ( bareMesh.points.numberOfColumns() ),
         numberBoundaryVertices ( 0 ),
         numberPoints           ( bareMesh.points.numberOfColumns() ),
         numberBoundaryPoints   ( 0 ),
         numberEdges            ( bareMesh.edges.numberOfColumns() ),
         numberBoundaryEdges    ( 0 ),
         numberFaces            ( bareMesh.numBoundaryFaces ),
         numberBoundaryFaces    ( 0 ),
         numberElements         ( bareMesh.elements.numberOfColumns() );
    UInt numberStoredFaces      ( bareMesh.faces.numberOfColumns() );

    std::vector<FaceHelp> faceHelp;
    typename std::vector<FaceHelp>::iterator faceHelpIterator;

    std::stringstream discardedLog;
    std::ostream& oStr = verbose ? std::cout : discardedLog;

    if ( verbose ) std::cout << "Converting bare mesh" << std::endl;

    // Be a little verbose
    switch ( GeoShape::S_shape )
    {
        case HEXA:
            if ( GeoShape::S_numPoints > 8 ) std::cout << "Quadratic Hexa mesh" << std::endl;
            else                             std::cout << "Linear Hexa mesh" << std::endl;
            break;

        case TETRA:
            if ( GeoShape::S_numPoints > 4 ) std::cout << "Quadratic Tetra mesh" << std::endl;
            else                             std::cout << "Linear Tetra Mesh" << std::endl;
            break;

        default:
            ERROR_MSG( "Current version of convertBareMesh only accepts TETRA and HEXA" );
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
         << "Number of Volumes         = "  << std::setw( 10 ) << numberElements         << std::endl;

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

    mesh.setMaxNumVolumes      ( numberElements, true );
    mesh.setMaxNumGlobalVolumes( numberElements );

    mesh.setMarkerID           ( bareMesh.regionMarkerID );

    typedef typename RegionMesh<GeoShape, MC>::point_Type  point_Type;
    typedef typename RegionMesh<GeoShape, MC>::edge_Type   edge_Type;
    typedef typename RegionMesh<GeoShape, MC>::face_Type   face_Type;
    typedef typename RegionMesh<GeoShape, MC>::volume_Type volume_Type;

    point_Type  * pointerPoint;
    edge_Type   * pointerEdge;
    face_Type   * pointerFace;
    volume_Type * pointerVolume;
    // addPoint()/Face()/Edge() returns a reference to the last stored point
    // I use that information to set all point info, by using a pointer.

    UInt count ( 0 );
    Int  ibc;

    // To account for internal faces
    if ( numberStoredFaces > numberBoundaryFaces )
    {
        faceHelp.resize( numberStoredFaces - numberBoundaryFaces );
        faceHelpIterator = faceHelp.begin();

        oStr << "WARNING: The mesh file (apparently) contains "
             << numberStoredFaces - numberBoundaryFaces << " internal faces" << std::endl;

    }

    // reading vertices
    for ( i = 0; i < numberVertices; i++ )
    {
        x = bareMesh.points( 0, i );
        y = bareMesh.points( 1, i );
        z = bareMesh.points( 2, i );
        ibc = bareMesh.pointsMarkers[ i ];

        bool isOnBoundary ( iSelect(markerID_Type(ibc)) );
        if ( isOnBoundary )
        {
            ++count;
        // Boundary point. Boundary switch set by the mesh method.
        }

        pointerPoint = &mesh.addPoint( isOnBoundary, true );

        pointerPoint->setId( i );
        pointerPoint->x() = x;
        pointerPoint->y() = y;
        pointerPoint->z() = z;
        pointerPoint->setMarkerID( markerID_Type( ibc ) );
        pointerPoint->setFlag( EntityFlags::VERTEX );
    }

    oStr << "Vertices read " << std::endl;
    oStr << "Size of the node storage is "
         << numberVertices * sizeof( point_Type ) / 1024. / 1024. << " MB" << std::endl;
    done++;

    if ( count != numberBoundaryVertices )
    {
        std::cerr << "Number boundary points inconsistent!" << std::endl;
    }

    oStr << "Reading boundary faces " << std::endl;

    // reading boundary faces
    for ( i = 0; i < numberStoredFaces; i++ )
    {
        for ( UInt j = 0; j < GeoShape::GeoBShape::S_numPoints; j++ )
        {
            p[ j ] = bareMesh.faces( j, i );
        }
        ibc = bareMesh.facesMarkers[ i ];

        bool isOnBoundary ( iSelect( markerID_Type( ibc ) ) );
        if ( isOnBoundary )
        {
            pointerFace = &( mesh.addFace( isOnBoundary ) ); // only boundary faces
            for ( UInt j = 0; j < GeoShape::GeoBShape::S_numPoints; j++ )
            {
                pointerFace->setPoint( j, mesh.point( p[ 0 ] ) );
            }
            pointerFace->setId( mesh.faceList.size() - 1 );
            pointerFace->setMarkerID( markerID_Type( ibc ) );
        }
        else // additional faces
        {
            for ( UInt j = 0; j < GeoShape::GeoBShape::S_numPoints; j++ )
            {
                faceHelpIterator->i[ j ] = p[ j ];
            }
            faceHelpIterator->ibc = ibc;
            ++faceHelpIterator;
        }
    }

    for ( faceHelpIterator = faceHelp.begin();
          faceHelpIterator != faceHelp.end(); ++faceHelpIterator )
    {
        for ( UInt j = 0; j < GeoShape::GeoBShape::S_numPoints; j++ )
        {
            p[ j ]  = faceHelpIterator->i[ j ];
        }
        ibc = faceHelpIterator->ibc;
        pointerFace  = &( mesh.addFace( false ) ); // INTERNAL FACE
        pointerFace->setId( mesh.faceList.size() - 1 );
        pointerFace->setMarkerID( markerID_Type( ibc ) );
        for ( UInt j = 0; j < GeoShape::GeoBShape::S_numPoints; j++ )
        {
            pointerFace->setPoint( j, mesh.point( p[ j ] ) ); // set face conn.
        }
    }

    oStr << "Boundary faces read " << std::endl;
    done++;


    oStr << "Reading boundary edges " << std::endl;

    for ( i = 0; i < numberBoundaryEdges; i++ )
    {
        for ( UInt j = 0; j < GeoShape::GeoBShape::GeoBShape::S_numPoints; j++ )
            p[ j ] = bareMesh.edges( j, i );
        ibc = bareMesh.edgesMarkers[ i ];
        pointerEdge = &mesh.addEdge( true ); // Only boundary edges.
        pointerEdge->setId( i );
        pointerEdge->setMarkerID( markerID_Type( ibc ) );
        for ( UInt j = 0; j < GeoShape::GeoBShape::GeoBShape::S_numPoints; j++ )
            pointerEdge->setPoint( j, mesh.point( p[ j ] ) ); // set edge conn.
    }
    oStr << "Boundary edges read " << std::endl;
    done++;

    count = 0;
    oStr << "Reading volumes " << std::endl;

    for ( i = 0; i < numberElements; i++ )
    {
        for ( UInt j = 0; j < GeoShape::S_numPoints; j++ )
            p[ j ] = bareMesh.elements( j, i );
        ibc = bareMesh.elementsMarkers[ i ];

        pointerVolume = &mesh.addVolume();
        pointerVolume->setId( i );
        for ( UInt j = 0; j < GeoShape::S_numPoints; j++ )
            pointerVolume->setPoint( j, mesh.point( p[ j ] ) );
        pointerVolume->setMarkerID( markerID_Type( ibc ) );
        count++;
    }
    oStr << "size of the volume storage is " << sizeof( volume_Type ) * count / 1024. / 1024.
         << " MB" << std::endl;
    oStr << count << " Volume elements read" << std::endl;
    done++;

    // Test the mesh
    Switch sw;

    // this if is the verbose version
    if ( !checkMesh3D( mesh, sw, true, verbose, oStr, std::cerr, oStr ) )
    {
        std::abort();
    }

    // This part is to build a P2 mesh from a P1 geometry
    if ( GeoShape::S_shape == TETRA && GeoShape::S_numPoints > 4 )
    {
        MeshUtility::p2MeshFromP1Data( mesh );
    }

    Real vols[ 3 ];
    getVolumeFromFaces( mesh, vols, oStr );

    oStr << "Volume enclosed by the mesh computed by integration on boundary faces"  << std::endl;
    oStr << "INT(X)     INT(Y)      INT(Z)     <- they should be equal and equal to" << std::endl
         << "                                     the voulume enclosed by the mesh"  << std::endl;
    oStr << vols[ 0 ] << "     " << vols[ 1 ] << "     " << vols[ 2 ] << std::endl;
    oStr << "Boundary faces are defining a closed surface if "
         << testClosedDomain( mesh, oStr ) << std::endl
         << "is (almost) zero" << std::endl;

    bareMesh.clear();

    return done == 4 ;
}// Function convertBareMesh


//! readINRIAMeshFile - reads .mesh meshes.
/*!
  @param bareMesh, the bareMesh data structure to fill in.
  @param fileName, the name of the mesh file  to read.
  @param regionFlag, the identifier for the region.
  @param verbose, setting it as true, the output is verbose (the default is false).
  @param iSelect,
  @return true if everything went fine, false otherwise.
*/
template <typename GeoShape>
bool
readINRIAMeshFile( RegionMeshBare<GeoShape> & bareMesh,
                   std::string const &        fileName,
                   markerID_Type              regionFlag,
                   bool                       verbose = false,
                   InternalEntitySelector     iSelect = InternalEntitySelector() )
{
    const int idOffset = 1; //IDs in INRIA meshes start from 1

    std::string line, faceName, volumeName;

    Real x, y, z;
    ID buffer;

    UInt done = 0;
    UInt numberVertices( 0 ), numberBoundaryVertices( 0 ),
         numberFaces   ( 0 ), numberBoundaryFaces   ( 0 ),
         numberPoints  ( 0 ), numberBoundaryPoints  ( 0 ),
         numberEdges   ( 0 ), numberBoundaryEdges   ( 0 );
    UInt numberVolumes ( 0 );
    UInt numberStoredFaces( 0 );

    std::stringstream discardedLog;

    ReferenceShapes shape( NONE );

    std::ostream& oStr = verbose ? std::cout : discardedLog;

    // open stream to read header

    std::ifstream hstream( fileName.c_str() );

    if ( verbose )
    {
        std::cout << "Reading from file " << fileName << std::endl;
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
        if ( verbose ) std::cout << "Linear Hexa mesh" << std::endl;
        numberPoints         =  numberVertices;
        numberBoundaryPoints = numberBoundaryVertices;
        faceName = "Quadrilaterals";
        volumeName = "Hexahedra";
        break;

    case TETRA:
        if ( GeoShape::S_numPoints == 6 )
        {
            if ( verbose ) std::cout << "Quadratic Tetra mesh (from linear geometry)" << std::endl;
            numberPoints         = numberVertices + numberEdges;
            numberBoundaryPoints = numberBoundaryVertices + numberBoundaryEdges;
        }
        else if ( GeoShape::S_numPoints == 4 )
        {
            if ( verbose ) std::cout << "Linear Tetra Mesh" << std::endl;

            numberPoints         = numberVertices;
            numberBoundaryPoints = numberBoundaryVertices;
            numberBoundaryEdges  = ( Int ( numberBoundaryVertices + numberBoundaryFaces - Int( 2 ) ) > 0 ?
                                         ( numberBoundaryVertices + numberBoundaryFaces - 2 ) : 0 );
        }
        else
        {
            ASSERT( 0, "mesh type not supported" );
        }

        faceName = "Triangles";
        volumeName = "Tetrahedra";
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
    bareMesh.points.reshape   ( 3, numberPoints );
    bareMesh.pointsMarkers.resize ( numberPoints );

    bareMesh.edges.reshape    ( GeoShape::GeoBShape::GeoBShape::S_numPoints, numberBoundaryEdges );
    bareMesh.edgesMarkers.resize ( numberBoundaryEdges );

    bareMesh.numBoundaryFaces = numberBoundaryFaces;
    bareMesh.faces.reshape    ( GeoShape::GeoBShape::S_numPoints, numberStoredFaces );
    bareMesh.facesMarkers.resize ( numberStoredFaces );

    bareMesh.elements.reshape ( GeoShape::S_numPoints, numberVolumes );
    bareMesh.elementsMarkers.resize ( numberVolumes );

    bareMesh.regionMarkerID = regionFlag;

    UInt count = 0;
    Int  ibc;

    // To account for internal faces
    if ( numberStoredFaces > numberBoundaryFaces )
    {
        oStr << "WARNING: The mesh file (apparently) contains "
             << numberStoredFaces - numberBoundaryFaces << " internal faces" << std::endl;
    }

    while ( nextGoodLine( myStream, line ).good() )
    {
        if ( line.find( "Vertices" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), myStream );

            for ( UInt i = 0; i < numberVertices; i++ )
            {
                myStream >> x >> y >> z >> ibc;

                if ( !iSelect(markerID_Type(ibc)))
                {
                    ++count;
                }

                bareMesh.points( 0, i ) = x;
                bareMesh.points( 1, i ) = y;
                bareMesh.points( 2, i ) = z;
                bareMesh.pointsMarkers[ i ] = ibc;
            }
            done++;

            if ( count != numberBoundaryVertices )
              {
                std::cerr << "Number boundary points inconsistent!" << std::endl;
              }
        }

        if ( line.find( faceName ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), myStream );
            oStr << "Reading boundary faces " << std::endl;

            for ( UInt i = 0; i < numberStoredFaces; i++ )
            {
                for( UInt k = 0; k < GeoShape::GeoBShape::S_numPoints; k++ )
                {
                    myStream >> buffer;
                    bareMesh.faces( k, i ) = buffer - idOffset;
                }
                bareMesh.facesMarkers[ i ] = ibc;
            }

            oStr << "Boundary faces read " << std::endl;
            done++;
        }

        if ( line.find( "Edges" ) != std::string::npos )
        {
            nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), myStream );
            oStr << "Reading boundary edges " << std::endl;

            for ( UInt i = 0; i < numberBoundaryEdges; i++ )
            {
                for( UInt k = 0; k < GeoShape::GeoBShape::GeoBShape::S_numPoints; k++ )
                {
                    myStream >> buffer;
                    bareMesh.edges( k, i ) = buffer - idOffset;
                }
                bareMesh.edgesMarkers[ i ] = ibc;
            }
            oStr << "Boundary edges read " << std::endl;
            done++;
        }

        if ( line.find( volumeName ) != std::string::npos )
        {
            count = 0;
            nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), myStream );
            oStr << "Reading volumes " << std::endl;

            for ( UInt i = 0; i < numberVolumes; i++ )
            {
                for( UInt k = 0; k < GeoShape::S_numPoints; k++ )
                {
                    myStream >> buffer;
                    bareMesh.elements( k, i ) = buffer - idOffset;
                }
                bareMesh.elementsMarkers[ i ] = ibc;
                count++;
            }
            oStr << count << " Volume elements read" << std::endl;
            done++;
        }
    }

    myStream.close();
    return done == 4 ;
}// Function readINRIAMeshFile

} // Namespace LifeV

#endif // _IMPORTERMESH3DNEW_HH_
