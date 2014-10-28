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
@brief Settings - File handling the solution of the FSI problem

@author Davide Forti <davide.forti@epfl.ch>
@date 28-10-2014

@maintainer Simone Deparis <simone.deparis@epfl.ch>
*/


#ifndef FSIHANDLER_H
#define FSIHANDLER_H

#include <lifev/core/LifeV.hpp>

// datafile
#include <lifev/core/filter/GetPot.hpp>

// includes for matrices and vector
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

// mesh and partitioner
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

namespace LifeV
{

//! FSIHandler - File handling the solution of the FSI problem
/*!
*  @author Davide Forti
*/

class FSIHandler
{
public:

    typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;
    
    typedef boost::shared_ptr<GetPot> datafilePtr_Type;
    
    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    
    typedef boost::shared_ptr<MeshData> meshDataPtr_Type;
    
    typedef boost::shared_ptr<MeshPartitioner<mesh_Type> > meshPartitionerPtr_Type;
    
    //! Constructor
    FSIHandler(const commPtr_Type& communicator);

    //! Destructor
    ~FSIHandler();
    
    void setDatafile(const GetPot& dataFile);

    void readMeshes( );
    
    void partitionMeshes ( );
    
//@}

private:

    commPtr_Type M_comm;
    
    GetPot M_datafile;
    
    meshDataPtr_Type M_meshDataFluid;
    meshPtr_Type M_fluidMesh;
    meshPtr_Type M_fluidLocalMesh;
    meshPartitionerPtr_Type M_fluidPartitioner;
    
    meshDataPtr_Type M_meshDataStructure;
    meshPtr_Type M_structureMesh;
    meshPtr_Type M_structureLocalMesh;
    meshPartitionerPtr_Type M_structurePartitioner;
    
};
    
} // end namespace LifeV

#endif // end SETTINGS_H
