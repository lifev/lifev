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
    @brief   File containing a class for handling spatial discretization.

    @author M.A. Fernandez
    @date 01-2003
    @author Cristiano Malossi<cristiano.malossi@epfl.ch>
    @date 06-2009
    @author Gilles Fourestey
    @date 06-2010
    @contributor Lucia Mirabella <lucia.mirabell@gmail.com>
    @maintainer Lucia Mirabella <lucia.mirabell@gmail.com>

 */

#include <lifeconfig.h>
#include <life/lifemesh/dataMesh.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
DataMesh::DataMesh( ):
        M_meshDir  ( "./" ),
        M_meshFile ( "mesh.mesh" ),
        M_meshType ( ".mesh" ),
        M_verbose   ( false )
{}

DataMesh::DataMesh( const GetPot& dataFile, const std::string& section ):
        M_meshDir  (),
        M_meshFile (),
        M_meshType (),
        M_verbose   ()
{
    setup( dataFile, section );
}

DataMesh::DataMesh( const DataMesh& dataMesh ):
        M_meshDir    ( dataMesh.M_meshDir ),
        M_meshFile   ( dataMesh.M_meshFile ),
        M_meshType   ( dataMesh.M_meshType ),
        M_verbose     ( dataMesh.M_verbose )
{}

// ===================================================
// Methods
// ===================================================
void
DataMesh::setup( const GetPot& dataFile, const std::string& section )
{
    M_meshDir  = dataFile( ( section + "/mesh_dir"  ).data(), "./" );
    M_meshFile = dataFile( ( section + "/mesh_file" ).data(), "mesh.mesh" );
    M_meshType = dataFile( ( section + "/mesh_type" ).data(), ".mesh" );
    M_verbose   = dataFile( ( section + "/verbose"   ).data(), false );
}

void DataMesh::showMe( std::ostream& output ) const
{
    output << "\n*** DataMesh: values for user-defined data\n\n";

    output << "mesh_dir   = " << M_meshDir  << std::endl;
    output << "mesh_file  = " << M_meshFile << std::endl;
    output << "mesh_type  = " << M_meshType << std::endl;
}

}
