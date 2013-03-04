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

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/MeshData.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MeshData::MeshData( ) :
    M_meshDir  ( "./" ),
    M_meshFile ( "mesh.mesh" ),
    M_meshType ( "structured" ),
    M_order ( "P1" ),
    M_verbose   ( false )
{}

MeshData::MeshData ( const GetPot& dataFile, const std::string& section ) :
    M_meshDir  (),
    M_meshFile (),
    M_meshType (),
    M_order    (),
    M_verbose   ()
{
    setup ( dataFile, section );
}

MeshData::MeshData ( const MeshData& meshData ) :
    M_meshDir    ( meshData.M_meshDir ),
    M_meshFile   ( meshData.M_meshFile ),
    M_meshType   ( meshData.M_meshType ),
    M_order      ( meshData.M_order ),
    M_verbose     ( meshData.M_verbose )
{}

// ===================================================
// Methods
// ===================================================
void
MeshData::setup ( const GetPot& dataFile, const std::string& section )
{
    M_meshDir  = dataFile ( ( section + "/mesh_dir"  ).data(), "./" );
    M_meshFile = dataFile ( ( section + "/mesh_file" ).data(), "mesh.mesh" );
    M_meshType = dataFile ( ( section + "/mesh_type" ).data(), "structured" );
    M_order    = dataFile ( ( section + "/mesh_order"   ).data(), "P1" );
    M_verbose   = dataFile ( ( section + "/verbose"   ).data(), false );
}

void
MeshData::setup ( const Teuchos::ParameterList& meshParameters )
{
    Teuchos::ParameterList defaultParameters;
    defaultParameters.set ("mesh_dir", "./");
    defaultParameters.set ("mesh_file", "mesh.mesh");
    defaultParameters.set ("mesh_type", ".mesh");
    defaultParameters.set ("order", "P1");
    defaultParameters.set ("verbose", false);
    defaultParameters.setParameters (meshParameters);

    M_meshDir = defaultParameters.get<std::string> ("mesh_dir");
    M_meshFile = defaultParameters.get<std::string> ("mesh_file");
    M_meshType = defaultParameters.get<std::string> ("mesh_type");
    M_order = defaultParameters.get<std::string> ("order");
    M_verbose = defaultParameters.get<bool> ("verbose");
}

void MeshData::showMe ( std::ostream& output ) const
{
    output << "\n*** MeshData: values for user-defined data\n\n";

    output << "mesh_dir   = " << M_meshDir  << std::endl;
    output << "mesh_file  = " << M_meshFile << std::endl;
    output << "mesh_type  = " << M_meshType << std::endl;
    output << "mesh_order  = " << M_order << std::endl;
}

}
