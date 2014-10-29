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

// solvers
#include <lifev/navier_stokes/solver/NavierStokesSolver.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialLinear.hpp>
/*
#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialNonLinear.hpp>
#include <lifev/structure/solver/isotropic/ExponentialMaterialNonLinear.hpp>
#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialNonLinearPenalized.hpp>
#include <lifev/structure/solver/isotropic/SecondOrderExponentialMaterialNonLinear.hpp>
#include <lifev/structure/solver/isotropic/NeoHookeanMaterialNonLinear.hpp>
*/
#include <lifev/fsi/solver/HarmonicExtensionSolver.hpp>

// Expression template FE space
#include <lifev/eta/fem/ETFESpace.hpp>

// time advance for the structure
#include <lifev/core/fem/TimeAdvance.hpp>
#include <lifev/core/fem/TimeAdvanceNewmark.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include <lifev/core/fem/TimeAndExtrapolationHandler.hpp>

namespace LifeV
{

//! FSIHandler - File handling the solution of the FSI problem
/*!
*  @author Davide Forti
*/

class FSIHandler
{
public:

	// Public typedefs

    typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;
    
    typedef boost::shared_ptr<GetPot> datafilePtr_Type;
    
    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    
    typedef boost::shared_ptr<MeshData> meshDataPtr_Type;
    
    typedef boost::shared_ptr<MeshPartitioner<mesh_Type> > meshPartitionerPtr_Type;
    
    typedef MapEpetra map_Type;
	typedef boost::shared_ptr<map_Type> mapPtr_Type;
    
    typedef FESpace< mesh_Type, map_Type > FESpace_Type;
    typedef boost::shared_ptr<FESpace_Type> FESpacePtr_Type;

    typedef ETFESpace< mesh_Type, map_Type, 3, 3 > solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type> solidETFESpacePtr_Type;

    typedef boost::shared_ptr<TimeAdvance<vector_Type> > timeAdvancePtr_Type;

    typedef BCHandler bc_Type;
    typedef boost::shared_ptr<BCHandler> bcPtr_Type;

	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;

    //! Constructor
    FSIHandler(const commPtr_Type& communicator);

    //! Destructor
    ~FSIHandler();
    
    void setDatafile(const GetPot& dataFile);

    void readMeshes( );
    
    void partitionMeshes ( );
    
    void setup ( );
    
    void setBoundaryConditions ( const bcPtr_Type& fluidBC, const bcPtr_Type& structureBC, const bcPtr_Type& aleBC);

    // update all the bc handlers
    void updateBoundaryConditions( );

    void initializeTimeAdvance ( );

//@}

private:

    void createStructureFESpaces ( );

    void createAleFESpace();

    // update the bc handler
    void updateBCHandler( bcPtr_Type & bc );

    //! communicator
    commPtr_Type M_comm;
    
    //! datafile
    GetPot M_datafile;
    
    // members for the fluid mesh
    meshDataPtr_Type M_meshDataFluid;
    meshPtr_Type M_fluidMesh;
    meshPtr_Type M_fluidLocalMesh;
    meshPartitionerPtr_Type M_fluidPartitioner;
    
    // members for the structure mesh
    meshDataPtr_Type M_meshDataStructure;
    meshPtr_Type M_structureMesh;
    meshPtr_Type M_structureLocalMesh;
    meshPartitionerPtr_Type M_structurePartitioner;
    
    // members for the fluid, the structura and the ALE fineite element spaces
    FESpacePtr_Type M_velocityFESpace;
    FESpacePtr_Type M_pressureFESpace;
    FESpacePtr_Type M_displacementFESpace;
    FESpacePtr_Type M_aleFESpace;

	solidETFESpacePtr_Type M_displacementETFESpace;
    
    // navier-stokes solver
    boost::shared_ptr<NavierStokesSolver> M_fluid;
    boost::shared_ptr<StructuralOperator<mesh_Type> > M_structure;
    boost::shared_ptr<StructuralConstitutiveLawData> M_dataStructure;
    boost::shared_ptr<HarmonicExtensionSolver<mesh_Type> > M_ale;
    
    // time advance for the structure
    timeAdvancePtr_Type M_structureTimeAdvance;
    boost::shared_ptr<TimeAndExtrapolationHandler> M_fluidTimeAdvance;

    // boundary conditions
    bcPtr_Type M_fluidBC;
    bcPtr_Type M_structureBC;
    bcPtr_Type M_aleBC;

	//! Displayer to print in parallel (only PID 0 will print)
	Displayer M_displayer;

	Real M_dt, M_t_zero, M_t_end;
};
    
} // end namespace LifeV

#endif // end SETTINGS_H
