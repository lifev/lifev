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
/*!
  \file dataMesh.h
  \author M.A. Fernandez
  \date 01/2003
  \version 1.0
 
  \brief File containing a class for handling spatial discretization with GetPot
 
*/
#ifndef _DATAMESH_H_
#define _DATAMESH_H_
#include <string>
#include <iostream>
#include "GetPot.hpp"
#include "lifeV.hpp"
#include "regionMesh3D.hpp"
#include "readMesh3D.hpp"

namespace LifeV
{

/*!
  \class DataMesh
 
  Base class which holds data concerning spatial discretization
 
*/
template <typename Mesh>
class DataMesh
{
public:

    //! Constructor
    /*!
      \param section the section in the data file
    */
    DataMesh( const GetPot& dfile, const std::string& section = "discretization" );

    //! Output
    virtual void showMe( std::ostream& c = std::cout ) const;

    //! The mesh
    Mesh& mesh();
    //! Virtual destructor
    virtual ~DataMesh();

protected:

    //! mesh
    std::string _mesh_dir;   // mesh dir
    std::string _mesh_file;  // mesh files
    std::string _mesh_type;  // mesh fil
    std::string _mesh_faces; // update all mesh faces
    std::string _mesh_edges; // update all mesh edges
    Mesh _mesh;       // the mesh

};


//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh>
DataMesh<Mesh>::
DataMesh( const GetPot& dfile, const std::string& section )
{
    _mesh_dir = dfile( ( section + "/mesh_dir" ).data(), "./" );
    _mesh_file = dfile( ( section + "/mesh_file" ).data(), "mesh.mesh" );
    _mesh_type = dfile( ( section + "/mesh_type" ).data(), ".mesh" );
    _mesh_faces = dfile( ( section + "/mesh_faces" ).data(), "boundary" );
    _mesh_edges = dfile( ( section + "/mesh_edges" ).data(), "boundary" );


    if ( _mesh_type == ".mesh" )
        readINRIAMeshFile( _mesh, _mesh_dir + _mesh_file, 1 ); // mesh readding
    else if ( _mesh_type == ".m++" )
        readMppFile( _mesh, _mesh_dir + _mesh_file, 1 );
    else
        ERROR_MSG( "Sorry, this mesh file can not be loaded" );

    if ( _mesh_edges == "all" )
        _mesh.updateElementEdges( true );
    else
        _mesh.updateElementEdges();
    if ( _mesh_faces == "all" )
        _mesh.updateElementFaces( true );
    else
        _mesh.updateElementFaces();

}

// Destructor
template <typename Mesh>
DataMesh<Mesh>::
~DataMesh()
{}


// Output
template <typename Mesh>
void DataMesh<Mesh>::
showMe( std::ostream& c ) const
{
    // mesh
    c << "mesh_dir   = " << _mesh_dir << std::endl;
    c << "mesh_file  = " << _mesh_file << std::endl;
    c << "mesh_type  = " << _mesh_type << std::endl;
    c << "mesh_edges = " << _mesh_edges << std::endl;
    c << "mesh_faces = " << _mesh_faces << std::endl;
}


// The mesh
template <typename Mesh>
Mesh& DataMesh<Mesh>::
mesh()
{
    return _mesh;
}
}
#endif
