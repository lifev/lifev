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
 @brief Settings - File handling the solution of the FSI problem
 
 @author Davide Forti <davide.forti@epfl.ch>
 @date 28-10-2014
 
@maintainer Simone Deparis <simone.deparis@epfl.ch>
*/

#include <lifev/core/LifeV.hpp>
#include <lifev/fsi_blocks/solver/FSIHandler.hpp>

namespace LifeV
{

FSIHandler::FSIHandler( const commPtr_Type& communicator ) :
M_comm ( communicator )
{
}

FSIHandler::~FSIHandler( )
{
}

void
FSIHandler::setDatafile( const GetPot& dataFile)
{
    M_datafile = dataFile;
}
    
void
FSIHandler::readMeshes( )
{
    M_fluidMesh.reset ( new mesh_Type ( M_comm ) );
    M_meshDataFluid.reset( new MeshData ( ) );
    M_meshDataFluid->setup (M_datafile, "fluid/space_discretization");
    readMesh (*M_fluidMesh, *M_meshDataFluid);
    
    M_structureMesh.reset ( new mesh_Type ( M_comm ) );
    M_meshDataStructure.reset( new MeshData ( ) );
    M_meshDataStructure->setup (M_datafile, "structure/space_discretization");
    readMesh (*M_structureMesh, *M_meshDataStructure);
}

void
FSIHandler::partitionMeshes( )
{
    M_fluidPartitioner.reset ( new MeshPartitioner< mesh_Type > (M_fluidMesh, M_comm) );
    M_fluidLocalMesh.reset ( new mesh_Type ( M_comm ) );
    M_fluidLocalMesh = M_fluidPartitioner->meshPartition();
    
    M_structurePartitioner.reset ( new MeshPartitioner< mesh_Type > (M_structureMesh, M_comm) );
    M_structureLocalMesh.reset ( new mesh_Type ( M_comm ) );
    M_structureLocalMesh = M_structurePartitioner->meshPartition();
}
    
void FSIHandler::setup ( )
{
    
}
    
    
}
