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

namespace
{

// local structure to store a bare info on a face
template <UInt n>
struct FacetStore
{
    ID points[ n ];
    markerID_Type bc;
};

}

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
template <typename GeoShapeType, typename MCType>
bool
convertBareMesh ( BareMesh<GeoShapeType> & bareMesh,
                  RegionMesh<GeoShapeType, MCType>&  mesh,
                  bool                       verbose = false,
                  InternalEntitySelector     iSelect = InternalEntitySelector() )
{
    typedef typename RegionMesh<GeoShapeType, MCType>::point_Type  point_Type;
    typedef typename RegionMesh<GeoShapeType, MCType>::ridge_Type   ridge_Type;
    typedef typename RegionMesh<GeoShapeType, MCType>::facet_Type   facet_Type;
    typedef typename RegionMesh<GeoShapeType, MCType>::element_Type element_Type;
    typedef std::vector<FacetStore<facet_Type::S_numPoints> > facetStore_Type;

    std::vector<UInt> points( element_Type::S_numVertices );

    UInt done = 0;
    UInt numberVertices         = bareMesh.points.numberOfColumns();
    UInt numberBoundaryVertices = bareMesh.numBoundaryPoints;
    UInt numberPoints           = bareMesh.points.numberOfColumns();
    UInt numberBoundaryPoints   = bareMesh.numBoundaryPoints;
    UInt numberStoredRidges     = bareMesh.ridges.numberOfColumns();
    UInt numberBoundaryRidges   = 0;
    UInt numberStoredFacets     = bareMesh.facets.numberOfColumns();
    UInt numberBoundaryFacets   = bareMesh.numBoundaryFacets;
    UInt numberElements         = bareMesh.elements.numberOfColumns();

    // Euler formulas
    UInt numberFacets = 2 * numberElements + ( numberBoundaryFacets / 2 );
    UInt numberRidges = ( 3 * numberBoundaryFacets - 2 * numberBoundaryVertices ) / 4 + numberVertices + numberElements;

    facetStore_Type facetHelp;
    typename facetStore_Type::iterator facetHelpIterator;

    std::stringstream discardedLog;
    std::ostream& oStr = verbose ? std::cout : discardedLog;

    if ( verbose ) std::cout << "Converting bare mesh" << std::endl;

    // Be a little verbose
    switch ( GeoShapeType::S_shape )
    {
        case HEXA:
            if     ( GeoShapeType::S_numPoints == 27 ) std::cout << "Quadratic Hexa mesh" << std::endl;
            else if( GeoShapeType::S_numPoints ==  8 ) std::cout << "Linear Hexa mesh" << std::endl;
            break;

        case TETRA:
            if     ( GeoShapeType::S_numPoints == 10 ) std::cout << "Quadratic Tetra mesh" << std::endl;
            else if( GeoShapeType::S_numPoints ==  4 ) std::cout << "Linear Tetra Mesh" << std::endl;
            break;

        default:
            ERROR_MSG( "Current version of convertBareMesh only accepts TETRA and HEXA" );
    }

    oStr << "Number of Vertices        = "  << std::setw( 10 ) << numberVertices         << std::endl
         << "Number of BVertices       = "  << std::setw( 10 ) << numberBoundaryVertices << std::endl
         << "Number of Faces           = "  << std::setw( 10 ) << numberFacets           << std::endl
         << "Number of Boundary Faces  = "  << std::setw( 10 ) << numberBoundaryFacets   << std::endl
         << "Number of Stored Faces    = "  << std::setw( 10 ) << numberStoredFacets     << std::endl
         << "Number of Edges           = "  << std::setw( 10 ) << numberRidges           << std::endl
         << "Number of Boundary Edges  = "  << std::setw( 10 ) << numberBoundaryRidges   << std::endl
         << "Number of Points          = "  << std::setw( 10 ) << numberPoints           << std::endl
         << "Number of Boundary Points = "  << std::setw( 10 ) << numberBoundaryPoints   << std::endl
         << "Number of Volumes         = "  << std::setw( 10 ) << numberElements         << std::endl;

    // Set all basic data structure

    // I store all Points
    mesh.setMaxNumPoints        ( numberPoints, true );
    mesh.setMaxNumGlobalPoints  ( numberPoints );
    mesh.setNumBPoints          ( numberBoundaryPoints );
    mesh.setNumVertices         ( numberVertices );
    mesh.setNumGlobalVertices   ( numberVertices );
    mesh.setNumBVertices        ( numberBoundaryVertices );
    // Only Boundary Edges (in a next version I will allow for different choices)
    mesh.setMaxNumRidges        ( numberBoundaryRidges );
    mesh.setMaxNumGlobalRidges  ( numberRidges );
    mesh.setNumRidges           ( numberRidges ); // Here the REAL number of edges (all of them)
    mesh.setNumBoundaryRidges   ( numberBoundaryRidges );
    // Only Boundary Faces
    mesh.setMaxNumFacets        ( numberStoredFacets );
    mesh.setMaxNumGlobalFacets  ( numberBoundaryFacets );
    mesh.setNumFacets           ( numberFacets ); // Here the REAL number of faces (all of them)
    mesh.setNumBoundaryFacets   ( numberBoundaryFacets );

    mesh.setMaxNumElements      ( numberElements, true );
    mesh.setMaxNumGlobalElements( numberElements );

    mesh.setMarkerID            ( bareMesh.regionMarkerID );

    // To account for internal faces
    if ( numberStoredFacets > numberBoundaryFacets )
    {
        facetHelp.resize( numberStoredFacets - numberBoundaryFacets );
        facetHelpIterator = facetHelp.begin();

        oStr << "WARNING: The mesh file (apparently) contains "
             << numberStoredFacets - numberBoundaryFacets << " internal faces" << std::endl;
    }

    // reading vertices
    UInt boundaryPointsCount = 0;
    for ( UInt i = 0; i < numberVertices; i++ )
    {
        Real const & x = bareMesh.points( 0, i );
        Real const & y = bareMesh.points( 1, i );
        Real const & z = bareMesh.points( 2, i );
        markerID_Type ibc = bareMesh.pointMarkers[ i ];

        bool isOnBoundary = !iSelect( ibc );
        if ( isOnBoundary )
        {
            boundaryPointsCount++;
        }

        point_Type* pointerPoint = &mesh.addPoint( isOnBoundary, true );
        pointerPoint->setId( i );
        pointerPoint->x() = x;
        pointerPoint->y() = y;
        pointerPoint->z() = z;
        pointerPoint->setMarkerID( ibc );
    }

    oStr << "Vertices read " << std::endl;
    oStr << "Size of the node storage is "
         << numberVertices * sizeof( point_Type ) / 1024. / 1024. << " MB" << std::endl;
    done++;

    if ( boundaryPointsCount != numberBoundaryVertices )
    {
        std::cerr << "Number boundary points inconsistent!" << std::endl;
    }

    oStr << "Reading boundary faces " << std::endl;

    // reading boundary faces
    for ( UInt i = 0; i < numberStoredFacets; i++ )
    {
        for ( UInt j = 0; j < facet_Type::S_numPoints; j++ )
        {
            points[ j ] = bareMesh.facets( j, i );
        }
        markerID_Type ibc = bareMesh.facetMarkers[ i ];

        bool isOnBoundary = !iSelect( ibc );
        if ( isOnBoundary )
        {
            facet_Type* pointerFace = &( mesh.addFacet( isOnBoundary ) ); // only boundary faces
            for ( UInt j = 0; j < facet_Type::S_numPoints; j++ )
            {
                pointerFace->setPoint( j, mesh.point( points[ j ] ) );
            }
            pointerFace->setId( mesh.faceList.size() - 1 );
            pointerFace->setMarkerID( ibc );
        }
        else // additional faces
        {
            for ( UInt j = 0; j < facet_Type::S_numPoints; j++ )
            {
                facetHelpIterator->points[ j ] = points[ j ];
            }
            facetHelpIterator->bc = ibc;
            ++facetHelpIterator;
        }
    }

    for ( facetHelpIterator = facetHelp.begin();
          facetHelpIterator != facetHelp.end(); ++facetHelpIterator )
    {
        for ( UInt j = 0; j < facet_Type::S_numPoints; j++ )
        {
            points[ j ]  = facetHelpIterator->points[ j ];
        }
        markerID_Type ibc = facetHelpIterator->bc;

        facet_Type* FacetP = &( mesh.addFacet( false ) ); // INTERNAL FACET
        FacetP->setId( mesh.faceList.size() - 1 );
        FacetP->setMarkerID( ibc );
        for ( UInt j = 0; j < facet_Type::S_numPoints; j++ )
        {
            FacetP->setPoint( j, mesh.point( points[ j ] ) ); // set face conn.
        }
    }

    oStr << "Boundary faces read " << std::endl;
    done++;


    oStr << "Reading boundary edges " << std::endl;

    for ( UInt i = 0; i < numberStoredRidges; i++ )
    {
        for ( UInt j = 0; j < ridge_Type::S_numPoints; j++ )
            points[ j ] = bareMesh.ridges( j, i );
        markerID_Type ibc = bareMesh.ridgeMarkers[ i ];
        ridge_Type* RidgeP = &mesh.addRidge( true ); // Only boundary edges.
        RidgeP->setId( i );
        RidgeP->setMarkerID( markerID_Type( ibc ) );
        for ( UInt j = 0; j < ridge_Type::S_numPoints; j++ )
            RidgeP->setPoint( j, mesh.point( points[ j ] ) ); // set edge conn.
    }
    oStr << "Boundary edges read " << std::endl;
    done++;

    oStr << "Reading volumes " << std::endl;

    for ( UInt i = 0; i < numberElements; i++ )
    {
        for ( UInt j = 0; j < element_Type::S_numPoints; j++ )
            points[ j ] = bareMesh.elements( j, i );
        markerID_Type ibc = bareMesh.elementMarkers[ i ];

        element_Type* ElementP = &mesh.addElement();
        ElementP->setId( i );
        for ( UInt j = 0; j < element_Type::S_numPoints; j++ )
            ElementP->setPoint( j, mesh.point( points[ j ] ) );
        ElementP->setMarkerID( ibc );
    }
    oStr << "size of the volume storage is " << sizeof( element_Type ) * numberElements / 1024. / 1024.
         << " MB" << std::endl;
    oStr << numberElements << " Volume elements read" << std::endl;
    done++;

    // Test the mesh
    Switch sw;

    // this if is the verbose version
    if ( !checkMesh3D( mesh, sw, true, verbose, oStr, std::cerr, oStr ) )
    {
        std::abort();
    }

    // This part is to build a P2 mesh from a P1 geometry
    if ( GeoShapeType::S_shape == TETRA && GeoShapeType::S_numPoints == 10 )
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
readINRIAMeshFile( BareMesh<GeoShape> & bareMesh,
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
    UInt numberPoints, numberBoundaryPoints;
    UInt numberVertices, numberBoundaryVertices;
    UInt numberEdges, numberBoundaryEdges;
    UInt numberBoundaryFaces, numberStoredFaces;
    UInt numberVolumes;

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

    // Set all basic data structure
    bareMesh.numBoundaryPoints = numberBoundaryPoints;
    bareMesh.points.reshape( 3, numberPoints );
    bareMesh.pointMarkers.resize( numberPoints );

    bareMesh.ridges.reshape    ( GeoShape::GeoBShape::GeoBShape::S_numPoints, numberBoundaryEdges );
    bareMesh.ridgeMarkers.resize( numberBoundaryEdges );

    bareMesh.numBoundaryFacets = numberBoundaryFaces;
    bareMesh.facets.reshape( GeoShape::GeoBShape::S_numPoints, numberStoredFaces );
    bareMesh.facetMarkers.resize( numberStoredFaces );

    bareMesh.elements.reshape( GeoShape::S_numPoints, numberVolumes );
    bareMesh.elementMarkers.resize( numberVolumes );

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
                bareMesh.pointMarkers[ i ] = ibc;
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
                    bareMesh.facets( k, i ) = buffer - idOffset;
                }
                myStream >> ibc;
                bareMesh.facetMarkers[ i ] = ibc;
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
                    bareMesh.ridges( k, i ) = buffer - idOffset;
                }
                myStream >> ibc;
                bareMesh.ridgeMarkers[ i ] = ibc;
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
                myStream >> ibc;
                bareMesh.elementMarkers[ i ] = ibc;
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
