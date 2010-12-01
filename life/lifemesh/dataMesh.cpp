//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2009-2010 EPFL, Politecnico di Milano

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a class for handling spatial discretization.
 *
 *  @author M.A. Fernandez
 *  @date 01/2003
 *  @version 1.0
 *
 *  @version 1.2
 *  @date 06/2009
 *  @author Cristiano Malossi<cristiano.malossi@epfl.ch>
 *
 *  @version 1.3
 *  @date 06/2010
 *  @author Gilles Fourestey
 */

#include <life/lifemesh/dataMesh.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
DataMesh::DataMesh( ):
        M_mesh_dir  ( "./" ),
        M_mesh_file ( "mesh.mesh" ),
        M_mesh_type ( ".mesh" ),
        M_verbose   ( false )
{}

DataMesh::DataMesh( const GetPot& dataFile, const std::string& section ):
        M_mesh_dir  (),
        M_mesh_file (),
        M_mesh_type (),
        M_verbose   ()
{
    setup( dataFile, section );
}

DataMesh::DataMesh( const DataMesh& dataMesh ):
        M_mesh_dir    ( dataMesh.M_mesh_dir ),
        M_mesh_file   ( dataMesh.M_mesh_file ),
        M_mesh_type   ( dataMesh.M_mesh_type ),
        M_verbose     ( dataMesh.M_verbose )
{}

// ===================================================
// Methods
// ===================================================
void
DataMesh::setup( const GetPot& dataFile, const std::string& section )
{
    M_mesh_dir  = dataFile( ( section + "/mesh_dir"  ).data(), "./" );
    M_mesh_file = dataFile( ( section + "/mesh_file" ).data(), "mesh.mesh" );
    M_mesh_type = dataFile( ( section + "/mesh_type" ).data(), ".mesh" );
    M_verbose   = dataFile( ( section + "/verbose"   ).data(), false );
}

void DataMesh::showMe( std::ostream& output ) const
{
    output << "\n*** DataMesh: values for user-defined data\n\n";

    output << "mesh_dir   = " << M_mesh_dir  << std::endl;
    output << "mesh_file  = " << M_mesh_file << std::endl;
    output << "mesh_type  = " << M_mesh_type << std::endl;
}

}
