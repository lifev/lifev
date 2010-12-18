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
 */

#include <life/lifefilters/readMesh3D.hpp>

namespace LifeV
{

// ===================================================
// Mpp mesh readers
// ===================================================

bool
readMppFileHead( std::ifstream & myStream,
                 UInt          & numberVertices,
                 UInt          & numberBoundaryVertices,
                 UInt          & numberBoundaryFaces,
                 UInt          & numberBoundaryEdges,
                 UInt          & numberVolumes )
{
    std::string line;

    Real x, y, z;

    Int ity, ity_id;

    UInt done = 0;
    UInt i, ibc;
    UInt p1, p2, p3;

    while ( next_good_line( myStream, line ).good() )
    {
        if ( line.find( "odes" ) != std::string::npos )
        {
            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );
            numberVertices = atoi( node_s );
            done++;
            numberBoundaryVertices = 0;

            for ( i = 0; i < numberVertices; i++ )
            {
                myStream >> x >> y >> z >> ity >> ibc;
                if ( ity != 3 )
                {

#ifndef OLDMPPFILE
                    myStream >> ibc;
#endif

                    numberBoundaryVertices++;
                }
            }
        }

        if ( line.find( "iangular" ) != std::string::npos )
        {
            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );
            numberBoundaryFaces = atoi( node_s );
            done++;

            for ( i = 0; i < numberBoundaryFaces; i++ )
            {
#ifdef OLDMPPFILE
                myStream >> p1 >> p2 >> p3 >> ity >> ibc;
#else

                myStream >> p1 >> p2 >> p3 >> ity >> ity_id >> ibc;
#endif
            }
        }

        if ( line.find( "oundary" ) != std::string::npos )
        {
            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );
            numberBoundaryEdges = atoi( node_s );

            for ( i = 0; i < numberBoundaryEdges; i++ )
            {
#ifdef OLDMPPFILE
                myStream >> p1 >> p2 >> ity >> ibc;
#else
                myStream >> p1 >> p2 >> ity >> ity_id >> ibc;
#endif
            }
            done++;
        }

        if ( line.find( "etrahedral" ) != std::string::npos )
        {
            std::string node_s = line.substr( line.find_last_of( ":" ) + 1 );
            numberVolumes = atoi( node_s );
            done++;
        }
    }

    return done == 4 ;
}// Function readMppFileHead

// ===================================================
// INRIA mesh readers
// ===================================================

Int 
nextIntINRIAMeshField( std::string const & line, 
                       std::istream      & myStream )
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
            return atoi( line );
          }
         }
    Int dummy;
    myStream >> dummy;
    
    return dummy;
}// Function nextIntINRIAMeshField

/*
 It Reads all basic info from INRIA MESH file
 so as to be able to properly dimension all arrays
*/
bool
readINRIAMeshFileHead( std::ifstream          & myStream,
                       UInt                   & numberVertices,
                       UInt                   & numberBoundaryVertices,
                       UInt                   & numberBoundaryFaces,
                       UInt                   & numberBoundaryEdges,
                       UInt                   & numberVolumes,
                       UInt                   & numberStoredFaces,
                       ReferenceShapes        & shape,
                       InternalEntitySelector   iSelect )
{
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

    while ( next_good_line( myStream, line ).good() )
    {
        if ( line.find( "MeshVersionFormatted" ) != std::string::npos )
        {
            idummy = nextIntINRIAMeshField( line.substr( line.find_last_of( "d" ) + 1 ), myStream );
            ASSERT_PRE0( idummy == 1, "I can read only formatted INRIA Mesh files, sorry" );
        }

        if ( line.find( "Dimension" ) != std::string:: npos )
        {
            idummy = nextIntINRIAMeshField( line.substr( line.find_last_of( "n" ) + 1 ), myStream );
            ASSERT_PRE0( idummy == 3, "I can read only 3D INRIA Mesh files, sorry" );
        }

        // I assume that internal vertices have their Ref value set to 0 (not clear from medit manual)
        if ( line.find( "Vertices" ) != std::string::npos )
        {
            numberVertices = nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), myStream );
            done++;
            numberBoundaryVertices = 0;
            isboundary.resize( numberVertices,false );

            for ( i = 0; i < numberVertices; ++i )
            {
                myStream >> x >> y >> z >> ibc;
                if ( ! iSelect( EntityFlag( ibc ) ) )
                {
                    numberBoundaryVertices++;
                    isboundary[ i ] = true;
                }
            }
        }

        if ( line.find( "Triangles" ) != std::string::npos )
        {
            ASSERT_PRE0( shape != HEXA, " Cannot have triangular faces in an HEXA INRIA  MESH" );
            shape = TETRA;
            numReadFaces = nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), myStream );
            numberBoundaryFaces = 0;

            done++;

            for ( UInt k = 0; k < numReadFaces; k++ )
            {
                myStream >> p1 >> p2 >> p3 >> ibc;
                if ( isboundary[ p1 - 1 ] && isboundary [ p2 - 1 ] && isboundary[ p3 - 1 ])
                {
                    if ( iSelect( EntityFlag( ibc ) ) )
                    {
                        std::cerr << "ATTENTION: Face "
                                  << p1 << " "
                                  << p2 << " "
                                  << p3 << " has all vertices on the boundary yet is marked as interior: "
                                  << ibc << std::endl;
                    }
                    ++numberBoundaryFaces;
                }
                else
                {
                    if ( !iSelect( EntityFlag( ibc ) ) )
                    {
                        std::cerr << "ATTENTION: Face "
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


        if ( line.find( "Quadrilaterals" ) != std::string::npos )
        {
            ASSERT_PRE0( shape != TETRA, " Cannot have quad faces in an TETRA INRIA MESH" );
            shape = HEXA;
            numReadFaces = nextIntINRIAMeshField( line.substr( line.find_last_of( "s" ) + 1 ), myStream );
            done++;
            numberBoundaryFaces = 0;

            for ( UInt k = 0; k < numReadFaces; k++ )
            {
               myStream >> p1 >> p2 >> p3 >> p4 >> ibc;
               if ( isboundary[ p1 - 1 ] && isboundary[ p2 - 1 ]
                    && isboundary[ p3 - 1 ] && isboundary[ p4 - 1 ] )
                {
                    if ( iSelect( EntityFlag( ibc ) ) )
                    {
                        std::cerr << "ATTENTION: Face "
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
        if ( line.find( "Tetrahedra" ) != std::string::npos )
        {
            ASSERT_PRE0( shape != HEXA, " Cannot have tetras  in a HEXA INRIA MESH" );
            shape = TETRA;
            numberVolumes = nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), myStream );
            done++;

            for ( i = 0; i < numberVolumes; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> p4 >> ibc;
            }
        }

        if ( line.find( "Hexahedra" ) != std::string::npos )
        {
            ASSERT_PRE0( shape != TETRA, " Cannot have Hexahedra in a TETRA INRIA MESH" );
            shape = HEXA;
            numberVolumes = nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), myStream );
            done++;

            for ( i = 0; i < numberVolumes; i++ )
            {
                myStream >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> p8 >> ibc;
            }
        }
        // I assume we are storing only boundary edges
        if ( line.find( "Edges" ) != std::string::npos )
        {
            numberBoundaryEdges = nextIntINRIAMeshField( line.substr( line.find_last_of( "a" ) + 1 ), myStream );
            done++;
            for ( i = 0; i < numberBoundaryEdges; i++ )
            {
                myStream >> p1 >> p2 >> ibc;
            }
        }
    }

    return true ;
}// Function readINRIAMeshFileHead

} // Namespace LifeV
