/********************************************************************************
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
********************************************************************************/

/**
 * @file   INRIAMeshParser.hpp
 * @brief  INRIA (Freefem++) mesh files (.mesh) reader.
 * @author Antonio Cervone <ant.cervone@gmail.com>
 * @date   10/2012
**/

#ifndef PARSER_INRIA_MESH_HPP__
#define PARSER_INRIA_MESH_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/BareMesh.hpp>
#include <lifev/core/mesh/InternalEntitySelector.hpp>

#include <fstream>

namespace LifeV
{

namespace MeshIO
{

namespace
{

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
nextIntINRIAMeshField ( std::string const& line,
                        std::istream&       myStream )
{
    /*
     first control if line has something.
     If so use atoi (the version from util_string.h)
     to extract the integer.
     Otherwise get if from the input stream
     */

    for ( std::string::const_iterator is = line.begin(); is != line.end(); ++is )
    {
        if ( *is != ' ' )
        {
            return atoi ( line );
        }
    }
    Int dummy;
    myStream >> dummy;

    return dummy;
}// Function nextIntINRIAMeshField

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
readINRIAMeshFileHead ( std::ifstream&           myStream,
                        UInt&                    numberVertices,
                        UInt&                    numberBoundaryVertices,
                        UInt&                    numberBoundaryFaces,
                        UInt&                    numberBoundaryEdges,
                        UInt&                    numberVolumes,
                        UInt&                    numberStoredFaces,
                        ReferenceShapes&         shape,
                        InternalEntitySelector   iSelect )
{
    const int idOffset = 1; //IDs in INRIA files start from 1

    numberVertices = 0;
    numberBoundaryVertices = 0;
    numberBoundaryFaces = 0;
    numberBoundaryEdges = 0;
    numberVolumes = 0;
    numberStoredFaces = 0;

    std::string line;

    Real x, y, z;

    Int idummy;

    UInt i, ibc;
    UInt done = 0;
    UInt p1, p2, p3, p4, p5, p6, p7, p8;
    UInt numReadFaces = 0;
    numberStoredFaces = 0;

    //shape = NONE;
    std::vector<bool> isboundary;
    //streampos start=myStream.tellg();

    while ( nextGoodLine ( myStream, line ).good() )
    {
        if ( line.find ( "MeshVersionFormatted" ) != std::string::npos )
        {
            idummy = nextIntINRIAMeshField ( line.substr ( line.find_last_of ( "d" ) + 1 ), myStream );
            ASSERT_PRE0 ( idummy == 1, "I can read only formatted INRIA Mesh files, sorry" );
        }

        if ( line.find ( "Dimension" ) != std::string:: npos )
        {
            idummy = nextIntINRIAMeshField ( line.substr ( line.find_last_of ( "n" ) + 1 ), myStream );
            ASSERT_PRE0 ( idummy == 3, "I can read only 3D INRIA Mesh files, sorry" );
        }

        // I assume that internal vertices have their Ref value set to 0 (not clear from medit manual)
        if ( line.find ( "Vertices" ) != std::string::npos )
        {
            numberVertices = nextIntINRIAMeshField ( line.substr ( line.find_last_of ( "s" ) + 1 ), myStream );
            done++;
            numberBoundaryVertices = 0;
            isboundary.resize ( numberVertices, false );

            for ( i = 0; i < numberVertices; ++i )
            {
                myStream >> x >> y >> z >> ibc;
                if ( ! iSelect ( markerID_Type ( ibc ) ) )
                {
                    numberBoundaryVertices++;
                    isboundary[ i ] = true;
                }
            }
        }

        if ( line.find ( "Triangles" ) != std::string::npos )
        {
            ASSERT_PRE0 ( shape != HEXA, " Cannot have triangular faces in an HEXA INRIA  MESH" );
            shape = TETRA;
            numReadFaces = nextIntINRIAMeshField ( line.substr ( line.find_last_of ( "s" ) + 1 ), myStream );
            numberBoundaryFaces = 0;

            done++;

            for ( UInt k = 0; k < numReadFaces; k++ )
            {
                myStream >> p1 >> p2 >> p3 >> ibc;
                if ( isboundary[ p1 - idOffset ] && isboundary [ p2 - idOffset ] && isboundary[ p3 - idOffset ])
                {
                    if ( iSelect ( markerID_Type ( ibc ) ) )
                    {
                        std::cerr << "ATTENTION: Face (1-based numbering) "
                                  << p1 << " "
                                  << p2 << " "
                                  << p3 << " has all vertices on the boundary yet is marked as interior: "
                                  << ibc << std::endl;
                    }
                    ++numberBoundaryFaces;
                }
                else
                {
                    if ( !iSelect ( markerID_Type ( ibc ) ) )
                    {
                        std::cerr << "ATTENTION: Face (1-based numbering) "
                                  << p1 << " "
                                  << p2 << " "
                                  << p3
                                  << " has vertices in the interior yet is marked as boundary: "
                                  << ibc << std::endl;
                    }
                }
            }
            numberStoredFaces = numReadFaces;
        }


        if ( line.find ( "Quadrilaterals" ) != std::string::npos )
        {
            ASSERT_PRE0 ( shape != TETRA, " Cannot have quad faces in an TETRA INRIA MESH" );
            shape = HEXA;
            numReadFaces = nextIntINRIAMeshField ( line.substr ( line.find_last_of ( "s" ) + 1 ), myStream );
            done++;
            numberBoundaryFaces = 0;

            for ( UInt k = 0; k < numReadFaces; k++ )
            {
                myStream >> p1 >> p2 >> p3 >> p4 >> ibc;
                if ( isboundary[ p1 - idOffset ] && isboundary[ p2 - idOffset ]
                        && isboundary[ p3 - idOffset ] && isboundary[ p4 - idOffset ] )
                {
                    if ( iSelect ( markerID_Type ( ibc ) ) )
                    {
                        std::cerr << "ATTENTION: Face (1-based numbering) "
                                  << p1 << " "
                                  << p2 << " "
                                  << p3 << " "
                                  << p4
                                  << " has all vertices on the boundary yet is marked as interior: "
                                  << ibc << std::endl;
                    }
                    ++numberBoundaryFaces;
                }

            }
            numberStoredFaces = numReadFaces;
        }
        // To cope with a mistake in INRIA Mesh files
        if ( line.find ( "Tetrahedra" ) != std::string::npos )
        {
            ASSERT_PRE0 ( shape != HEXA, " Cannot have tetras  in a HEXA INRIA MESH" );
            shape = TETRA;
            numberVolumes = nextIntINRIAMeshField ( line.substr ( line.find_last_of ( "a" ) + 1 ), myStream );
            done++;

            for ( i = 0; i < numberVolumes; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> p4 >> ibc;
            }
        }

        if ( line.find ( "Hexahedra" ) != std::string::npos )
        {
            ASSERT_PRE0 ( shape != TETRA, " Cannot have Hexahedra in a TETRA INRIA MESH" );
            shape = HEXA;
            numberVolumes = nextIntINRIAMeshField ( line.substr ( line.find_last_of ( "a" ) + 1 ), myStream );
            done++;

            for ( i = 0; i < numberVolumes; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> p8 >> ibc;
            }
        }
        // I assume we are storing only boundary edges
        if ( line.find ( "Edges" ) != std::string::npos )
        {
            numberBoundaryEdges = nextIntINRIAMeshField ( line.substr ( line.find_last_of ( "s" ) + 1 ), myStream );
            done++;
            for ( i = 0; i < numberBoundaryEdges; i++ )
            {
                myStream >> p1 >> p2 >> ibc;
            }
        }
    }

    return true ;
}// Function readINRIAMeshFileHead

}

//! INRIAMeshRead - reads .mesh meshes.
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
ReadINRIAMeshFile ( BareMesh<GeoShape>& bareMesh,
                    std::string const&         fileName,
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

    ReferenceShapes shape ( NONE );

    std::ostream& oStr = verbose ? std::cout : discardedLog;

    // open stream to read header

    std::ifstream hstream ( fileName.c_str() );

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

    if ( ! MeshIO::readINRIAMeshFileHead ( hstream, numberVertices, numberBoundaryVertices,
                                           numberBoundaryFaces, numberBoundaryEdges,
                                           numberVolumes, numberStoredFaces, shape, iSelect) )
    {
        std::cerr << " Error while reading INRIA mesh file headers" << std::endl;
        std::abort() ;
    }

    hstream.close();

    //Reopen the stream: I know it is stupid but this is how it goes
    std::ifstream myStream ( fileName.c_str() );

    if ( myStream.fail() )
    {
        std::cerr << " Error in readINRIAMeshFile = file " << fileName
                  << " not found or locked" << std::endl;
        std::abort();
    }

    ASSERT_PRE0 ( GeoShape::S_shape == shape, "INRIA Mesh file and mesh element shape is not consistent" );

    numberPoints         =  numberVertices;
    numberBoundaryPoints = numberBoundaryVertices;

    // Be a little verbose
    switch ( shape )
    {
        case HEXA:
            ASSERT_PRE0 ( GeoShape::S_numPoints == 8, "Sorry I can read only linear Hexa meshes" );
            if ( verbose )
            {
                std::cout << "Linear Hexa mesh" << std::endl;
            }
            faceName = "Quadrilaterals";
            volumeName = "Hexahedra";
            break;

        case TETRA:
            if ( GeoShape::S_numPoints == 4 )
            {
                if ( verbose )
                {
                    std::cout << "Linear Tetra Mesh" << std::endl;
                }
            }
            else if ( GeoShape::S_numPoints == 10 )
            {
                if ( verbose )
                {
                    std::cout << "Quadratic Tetra mesh (from linear geometry)" << std::endl;
                }
                numberPoints         += numberEdges;
                numberBoundaryPoints += numberBoundaryEdges;
            }
            else
            {
                ERROR_MSG ( "mesh type not supported" );
            }

            faceName = "Triangles";
            volumeName = "Tetrahedra";
            break;
        default:
            ERROR_MSG ( "Current version of INRIA Mesh file reader only accepts TETRA and HEXA" );
    }

    // Set all basic data structure
    bareMesh.numBoundaryPoints = numberBoundaryPoints;
    bareMesh.numVertices = numberVertices;
    bareMesh.numBoundaryVertices = numberBoundaryVertices;
    bareMesh.points.reshape ( 3, numberPoints );
    bareMesh.pointMarkers.resize ( numberPoints );
    bareMesh.pointIDs.resize ( numberPoints );

    bareMesh.ridges.reshape    ( GeoShape::GeoBShape::GeoBShape::S_numPoints, numberBoundaryEdges );
    bareMesh.ridgeMarkers.resize ( numberBoundaryEdges );
    bareMesh.ridgeIDs.resize ( numberBoundaryEdges );

    bareMesh.numBoundaryFacets = numberBoundaryFaces;
    bareMesh.facets.reshape ( GeoShape::GeoBShape::S_numPoints, numberStoredFaces );
    bareMesh.facetMarkers.resize ( numberStoredFaces );
    bareMesh.facetIDs.resize ( numberStoredFaces );

    bareMesh.elements.reshape ( GeoShape::S_numPoints, numberVolumes );
    bareMesh.elementMarkers.resize ( numberVolumes );
    bareMesh.elementIDs.resize ( numberVolumes );

    bareMesh.regionMarkerID = regionFlag;

    UInt count = 0;
    Int  ibc;

    // To account for internal faces
    if ( numberStoredFaces > numberBoundaryFaces )
    {
        oStr << "WARNING: The mesh file (apparently) contains "
             << numberStoredFaces - numberBoundaryFaces << " internal faces" << std::endl;
    }

    while ( nextGoodLine ( myStream, line ).good() )
    {
        if ( line.find ( "Vertices" ) != std::string::npos )
        {
            nextIntINRIAMeshField ( line.substr ( line.find_last_of ( "s" ) + 1 ), myStream );

            for ( UInt i = 0; i < numberVertices; i++ )
            {
                myStream >> x >> y >> z >> ibc;

                if ( !iSelect (markerID_Type (ibc) ) )
                {
                    ++count;
                }

                bareMesh.points ( 0, i ) = x;
                bareMesh.points ( 1, i ) = y;
                bareMesh.points ( 2, i ) = z;
                bareMesh.pointMarkers[ i ] = ibc;
                bareMesh.pointIDs[i] = i;
            }
            done++;

            if ( count != numberBoundaryVertices )
            {
                std::cerr << "Number boundary points inconsistent!" << std::endl;
            }
        }

        if ( line.find ( faceName ) != std::string::npos )
        {
            nextIntINRIAMeshField ( line.substr ( line.find_last_of ( "s" ) + 1 ), myStream );
            oStr << "Reading boundary faces " << std::endl;

            for ( UInt i = 0; i < numberStoredFaces; i++ )
            {
                for ( UInt k = 0; k < GeoShape::GeoBShape::S_numPoints; k++ )
                {
                    myStream >> buffer;
                    bareMesh.facets ( k, i ) = buffer - idOffset;
                }
                myStream >> ibc;
                bareMesh.facetMarkers[ i ] = ibc;
                bareMesh.facetIDs[ i ] = i;
            }

            oStr << "Boundary faces read " << std::endl;
            done++;
        }

        if ( line.find ( "Edges" ) != std::string::npos )
        {
            nextIntINRIAMeshField ( line.substr ( line.find_last_of ( "s" ) + 1 ), myStream );
            oStr << "Reading boundary edges " << std::endl;

            for ( UInt i = 0; i < numberBoundaryEdges; i++ )
            {
                for ( UInt k = 0; k < GeoShape::GeoBShape::GeoBShape::S_numPoints; k++ )
                {
                    myStream >> buffer;
                    bareMesh.ridges ( k, i ) = buffer - idOffset;
                }
                myStream >> ibc;
                bareMesh.ridgeMarkers[ i ] = ibc;
                bareMesh.ridgeIDs[ i ] = i;
            }
            oStr << "Boundary edges read " << std::endl;
            done++;
        }

        if ( line.find ( volumeName ) != std::string::npos )
        {
            count = 0;
            nextIntINRIAMeshField ( line.substr ( line.find_last_of ( "a" ) + 1 ), myStream );
            oStr << "Reading volumes " << std::endl;

            for ( UInt i = 0; i < numberVolumes; i++ )
            {
                for ( UInt k = 0; k < GeoShape::S_numPoints; k++ )
                {
                    myStream >> buffer;
                    bareMesh.elements ( k, i ) = buffer - idOffset;
                }
                myStream >> ibc;
                bareMesh.elementMarkers[ i ] = ibc;
                bareMesh.elementIDs[ i ] = i;
                count++;
            }
            oStr << count << " Volume elements read" << std::endl;
            done++;
        }
    }

    myStream.close();
    return done == 4 ;
}

} // GmshIO

} // LifeV

#endif // PARSER_INRIA_MESH_HPP__
