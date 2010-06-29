/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file dataMesh.cpp

*/

#include <life/lifemesh/dataMesh.hpp>

// ===================================================
// Constructors & Destructor
// ===================================================

namespace LifeV{


DataMesh::DataMesh( ):
        //M_mesh      ( new Mesh ),
        M_mesh_dir  ( "./" ),
        M_mesh_file ( "mesh.mesh" ),
        M_mesh_type ( ".mesh" ),
        M_verbose   ( false )
{}

DataMesh::DataMesh( const GetPot& dataFile, const std::string& section ):
        //M_mesh      ( new Mesh ),
        M_mesh_dir  ( dataFile( ( section + "/mesh_dir"  ).data(), "./" ) ),
        M_mesh_file ( dataFile( ( section + "/mesh_file" ).data(), "mesh.mesh" ) ),
        M_mesh_type ( dataFile( ( section + "/mesh_type" ).data(), ".mesh" ) ),
        M_verbose   ( dataFile( ( section + "/verbose" ).data(), false ) )
{
    //readMesh(dataFile);
}

DataMesh::DataMesh( const DataMesh& dataMesh ):
        //M_mesh        ( dataMesh.M_mesh ),
        M_mesh_dir    ( dataMesh.M_mesh_dir ),
        M_mesh_file   ( dataMesh.M_mesh_file ),
        M_mesh_type   ( dataMesh.M_mesh_type ),
        M_verbose     ( dataMesh.M_verbose )
{
    //M_mesh->updateElementEdges();
    //M_mesh->updateElementFaces();
}



// ===================================================
// Methods
// ===================================================
void
DataMesh::setup( const GetPot& dataFile, const std::string& section )
{
    M_mesh_dir  = dataFile( ( section + "/mesh_dir"  ).data(), "./" );
    M_mesh_file = dataFile( ( section + "/mesh_file" ).data(), "mesh.mesh" );
    M_mesh_type = dataFile( ( section + "/mesh_type" ).data(), ".mesh" );
    M_verbose   = dataFile( ( section + "/verbose" ).data(), 0 );

    //readMesh(dataFile);
}

void DataMesh::showMe( std::ostream& output ) const
{
    output << "mesh_dir   = " << M_mesh_dir << std::endl;
    output << "mesh_file  = " << M_mesh_file << std::endl;
    output << "mesh_type  = " << M_mesh_type << std::endl;
}

}
