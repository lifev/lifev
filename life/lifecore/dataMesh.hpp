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
  DataMesh(const GetPot& dfile, const string& section="discretization");
  
  //! Ouptut
  virtual void showMe(ostream& c=cout) const;

  //! The mesh
  Mesh& mesh();
  
  //! Virtual destructor
  virtual ~DataMesh();

 protected: 
 
  //! mesh 
  string _mesh_dir;  // mesh dir
  string _mesh_file; // mesh fil
  Mesh   _mesh;      // the mesh

};


//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh> 
DataMesh<Mesh>::
DataMesh(const GetPot& dfile, const string& section)
{
  _mesh_dir  = dfile((section+"/mesh_dir").data(),"./");
  _mesh_file = dfile((section+"/mesh_file").data(),"mesh.mesh");

#ifdef MESH_INRIA
  readINRIAMeshFile(_mesh, _mesh_dir+_mesh_file, 1); // mesh readding
#else // MESH_MOX
  readMppFile(_mesh,_mesh_dir+_mesh_file,1);
#endif
  _mesh.updateElementEdges();
  _mesh.updateElementFaces();
 
}

// Destructor
template <typename Mesh>
DataMesh<Mesh>::
~DataMesh() {}


// Output
template <typename Mesh>
void DataMesh<Mesh>::
showMe(ostream& c) const
{
  // mesh
  c << "mesh_dir  = " << _mesh_dir << endl; 
  c << "mesh_file = " << _mesh_file << endl; 
}


// The mesh
template <typename Mesh>
Mesh& DataMesh<Mesh>::
mesh() {
  return  _mesh;
}

#endif
