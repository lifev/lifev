/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Martin Prosi <martin.prosi@epfl.ch>
      Date: 2004-09-10

 Copyright (C) 2004 EPFL

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
/**
   \file ensight7Writer.hpp
   \author Martin Prosi <martin.prosi@epfl.ch>
   \date 2004-09-10
 */
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <string>

#include <life/lifecore/life.hpp>

namespace LifeV
{
class Coord
{
public:
    Coord()
    {
        init( 0 );
    };
    Coord &operator=( const Coord & A )
    {
        init( A.coord );
        return *this;
    };
    float &operator[] ( int i )
    {
        return coord[ i ];
    };

private:

    void init( const float* a )
    {
        for ( int i = 0;i < 3;++i )
            coord[ i ] = ( a != 0 ) ? a[ i ] : -100;
    };

    float coord[ 3 ];
};


bool fromensight7Mesh3D( UInt dimDof,
                         Vector & u,
                         Vector & p,
                         Real const & time,
                         std::string prefix = "./" )
{
	// store the values read from file
	std::vector<Coord> velocity( dimDof );
	std::vector<float> pressure( dimDof );

	char buffer[ 80 ], buf[ 80 ];

	// index of the file
	std::ostringstream findex;
	findex << std::setfill('0') << std::setw(3); 
	findex << ( time * 100 );

	// name of the file to be read
	std::string ifname( prefix + "velocity.res" );
	ifname += findex.str();

	std::fstream FileU( ifname.c_str(), std::ios::in | std::ios::binary );
	if (FileU.is_open())
	{
		std::cout << "\nreading file " << ifname << std::endl;
		// by construction, this is the header of the file
		strcpy( buffer, "velocity field  timestep 1" );
		FileU.read( ( char* ) & buf, sizeof( buffer ) );
		// read the vector of coordinates
		FileU.read( ( char* ) & velocity.front(), 
		            3 * velocity.size() * sizeof( float ) );

		FileU.close();		
	}
	else
	{
		std::cout << "Error opening file" << ifname << std::endl;
	}

  Coord nodevel;
//	std::cout << "\nprinting velocity values:" << std::endl;
	for( UInt i = 0; i < dimDof; ++i ) {
		nodevel = velocity[i];
		for( UInt d = 0; d < nDimensions; ++d ) {
//			std::cout << nodevel[d] << " " << std::flush;
			u[i + d*dimDof] = nodevel[d];
		}
	}
//	std::cout << std::endl;
	
	// name of the file to be read
	ifname = prefix + "pressure.res";
	ifname += findex.str();

	std::fstream FileP( ifname.c_str(), std::ios::in | std::ios::binary );
	if (FileP.is_open())
	{
		// by construction, this is the header of the file
		strcpy( buffer, "concentration distribution timestep " );
		FileP.read( ( char* ) & buf, sizeof( buffer ) );
		// read the vector of coordinates
		FileP.read( ( char* ) & pressure.front(), 
		            pressure.size() * sizeof( float ) );

		FileP.close();
	}
	else
	{
		std::cout << "Error opening file" << ifname << std::endl;
	}

	for( UInt i = 0; i < dimDof; ++i )
		p[i] = pressure[i];
	
	return 0;
}


template <typename RegionMesh3D>
bool outensight7Mesh3D( RegionMesh3D const & mesh,
                        PhysVectUnknown<Vector> const& u,
                        ScalUnknown<Vector> const & p,
                        Real const & time )
{

    std::vector<Coord> grid; // coordinates of Grid nodes
    grid.resize( mesh.numVertices() + mesh.numVolumes() );
    //  grid.resize(mesh.numVertices()+mesh.numEdges());

    for ( ID i = 0; i < mesh.numVertices(); i++ )
    {
        grid[ i ][ 0 ] = ( float ) mesh.point( i + 1 ).x();
        grid[ i ][ 1 ] = ( float ) mesh.point( i + 1 ).y();
        grid[ i ][ 2 ] = ( float ) mesh.point( i + 1 ).z();
    }
    //     std::cout << grid[i][0] << ", " << grid[i][1]  << ", " << grid[i][2] << std::endl; }

    // 6 additional mesh points for Tetra P2 (mittle points of edges)
    //  typename RegionMesh3D::EdgeType * pe=0;

    //  ID i1,i2;

    //  for(ID i=1; i <= mesh.numEdges(); i++){
    //     pe=& mesh.edge(i);
    //     i1=(pe->point(1)).id();
    //     i2=(pe->point(2)).id();
    //     grid[mesh.numVertices()+i-1][0]=(float)(mesh.point(i1).x()+mesh.point(i2).x())*0.5;
    //     grid[mesh.numVertices()+i-1][1]=(float)(mesh.point(i1).y()+mesh.point(i2).y())*0.5;
    //     grid[mesh.numVertices()+i-1][2]=(float)(mesh.point(i1).z()+mesh.point(i2).z())*0.5;}
    //     std::cout << grid[i][0] << ", " << grid[i][1]  << ", " << grid[i][2] << std::endl;}


    char buffer[ 80 ];
    std::vector<int> idnode, idelem;

    std::fstream File3( "test.geo", std::ios::out | std::ios::binary );

    strcpy( buffer, "C Binary" );
    File3.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "test cube - inria mesh" );
    File3.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "LinearTetra - elements" );
    File3.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "node id given" );
    File3.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "element id given" );
    File3.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "coordinates" );
    File3.write( ( char * ) & buffer, sizeof( buffer ) );
    int i = mesh.numVertices();
    //  int i=mesh.numVertices()+mesh.numVolumes();
    //  int i=mesh.numVertices() + mesh.numEdges();
    File3.write( ( char * ) & i, sizeof( int ) );
    for ( ID i = 1; i <= mesh.numVertices(); i++ )
    {  // load node id-vector for nodes
        //  for (ID i = 1; i <= mesh.numVertices()+mesh.numVolumes(); i++){  // load node id-vector for nodes
        //  for (ID i = 1; i <= mesh.numVertices()+mesh.numEdges(); i++){  // load node id-vector for nodes
        idnode.push_back( i );
    }
    File3.write( ( char * ) & idnode.front(), idnode.size() * sizeof( int ) );
    File3.write( ( char * ) & grid.front(), 3 * mesh.numVertices() * sizeof( float ) );
    //  File3.write((char *)&grid.front(),3*(mesh.numVertices()+mesh.numVolumes())*sizeof(float));
    //  File3.write((char *)&grid.front(),3*(mesh.numVertices()+mesh.numEdges())*sizeof(float));
    strcpy( buffer, "part 1" );
    File3.write( ( char * ) & buffer, sizeof( buffer ) );
    strcpy( buffer, "description line" );
    File3.write( ( char * ) & buffer, sizeof( buffer ) );
    if ( mesh.numLocalVertices() == 4 )
        strcpy( buffer, "tetra4" );
    //     strcpy(buffer,"tetra10");
    else
        strcpy( buffer, "hexa8" );
    File3.write( ( char * ) & buffer, sizeof( buffer ) );
    int e = mesh.storedVolumes();
    //  int e=4*mesh.storedVolumes();
    File3.write( ( char * ) & e, sizeof( int ) );
    for ( int i = 1; i <= e; i++ )
    {  // load node id-vector for elements
        idelem.push_back( i );
    }
    File3.write( ( char * ) & idelem.front(), idelem.size() * sizeof( int ) );

    // ------------------   Output P1 Elements (Tetra or Hexahedra)   -------------------------

    for ( ID k = 0; k < mesh.storedVolumes(); k++ )
    {
        for ( ID j = 0; j < mesh.numLocalVertices(); j++ )
        {
            int __id = mesh.volume( k + 1 ).point( j + 1 ).id();
            File3.write( ( char * ) & __id , sizeof( ID ) );
        }
    }

    // ------------------     Output P2 Elements (Tetra)     -----------------------------------

    // for( ID k = 0; k < mesh.storedVolumes(); k++){
    //      for (ID j = 0; j < mesh.numLocalVertices(); j++){
    // File3.write((char *)&mesh.volume(k+1).point(j+1).id(),sizeof(ID));
    //      }
    //      for (ID j = 0; j < mesh.numLocalEdges(); j++){
    //  ID i = mesh.localEdgeId(k+1,j+1)+mesh.numVertices();
    // File3.write((char *)&i,sizeof(ID));
    //      }
    //  }


    File3.close();

    std::cout << "output in ensight7 format" << std::endl;
    std::cout << "geometry file is test.geo" << std::endl;

    // read pressure von acsii-file ./Post/presQ1.bb and convert it into ensight7 binary format

    std::vector<float> pressure; // presure values at Grid nodes
    std::vector<Coord> velocity; // velocity values at Grid nodes

    std::ostringstream index;
    std::string name, vname;
    Coord nodevel;

    index << ( time * 100 );

    switch ( index.str().size() )
    {
    case 1:
        name = "pressure.res00" + index.str();
        vname = "velocity.res00" + index.str();
        break;
    case 2:
        name = "pressure.res0" + index.str();
        vname = "velocity.res0" + index.str();
        break;
    case 3:
        name = "pressure.res" + index.str();
        vname = "velocity.res" + index.str();
        break;
    }

    for ( ID i = 0; i < p.size(); i++ )
    {
        pressure.push_back( ( float ) ( p( i ) ) );
    }

    for ( ID i = 0; i < mesh.numVertices(); i++ )
    {
        //   for(ID i=0; i< mesh.numVertices()+mesh.numEdges(); i++){
        nodevel[ 0 ] = ( float ) u( i );
        nodevel[ 1 ] = ( float ) u( i + ( u.size() / 3 ) );
        nodevel[ 2 ] = ( float ) u( i + 2 * ( u.size() / 3 ) );
        velocity.push_back( nodevel );
    }

    std::fstream File4( name.c_str(), std::ios::out | std::ios::binary );

    strcpy( buffer, "concentration distribution timestep " );
    File4.write( ( char* ) & buffer, sizeof( buffer ) );

    File4.write( ( char* ) & pressure.front(), p.size() * sizeof( float ) );

    File4.close();

    std::fstream File5( vname.c_str(), std::ios::out | std::ios::binary );

    strcpy( buffer, "velocity field  timestep 1" );
    File5.write( ( char* ) & buffer, sizeof( buffer ) );

    File5.write( ( char* ) & velocity.front(), 3 * velocity.size() * sizeof( float ) );

    File5.close();

    std::cout << "result files are velocity.res*** and pressure.res*** " << std::endl;

    return true;

}

}
