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

 */

#ifndef NAVIERSTOKESSOLVERBLOCK_H
#define NAVIERSTOKESSOLVERBLOCK_H

#include <lifev/core/LifeV.hpp>

// includes for matrices and vector
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

// including data container for the fluid problem
#include <lifev/navier_stokes/solver/OseenData.hpp>

// includes for building the matrix graph
#include <Epetra_FECrsGraph.h>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/eta/expression/BuildGraph.hpp>

// classical FE space
#include <lifev/core/fem/FESpace.hpp>

// boundary conditions
#include <lifev/core/fem/BCManage.hpp>

// expression template finite element space
#include <lifev/eta/fem/ETFESpace.hpp>

// includes for the linear solver
#include <lifev/navier_stokes/solver/NavierStokesOperator.hpp>
#include <lifev/navier_stokes/solver/aSIMPLEOperator.hpp>

// utilities
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/util/Displayer.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

namespace LifeV
{

class NavierStokesSolver
{

public:

	typedef RegionMesh<LinearTetra> mesh_Type;
	typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

	typedef MapEpetra map_Type;
	typedef boost::shared_ptr<map_Type> mapPtr_Type;

	typedef MatrixEpetra<Real> matrix_Type;
	typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

	typedef VectorEpetra vector_Type;
	typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

	typedef GetPot dataFile_Type;
	typedef boost::shared_ptr<dataFile_Type> dataFilePtr_Type;

	typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;

	typedef Epetra_FECrsGraph graph_Type;
	typedef boost::shared_ptr<Epetra_FECrsGraph> graphPtr_Type;

	typedef ETFESpace<mesh_Type, map_Type, 3, 3 > ETFESpace_velocity;
	typedef ETFESpace<mesh_Type, map_Type, 3, 1 > ETFESpace_pressure;

	typedef boost::shared_ptr<BCHandler> bcPtr_Type;

	typedef Teuchos::ParameterList parameterList_Type;
	typedef boost::shared_ptr<parameterList_Type> parameterListPtr_Type;
    
    typedef boost::shared_ptr<Epetra_Operator> invOpPtr_Type;

	// Constructor
	NavierStokesSolver(const dataFile_Type dataFile, const commPtr_Type& communicator);

	// Destructor
	~NavierStokesSolver();

	// Setup
	void setup(const meshPtr_Type& mesh);

	// Assemble constant terms
	void buildSystem();

	// Update the convective term and the right hand side, sum contributions of block (0,0) in M_F
	void updateSystem( const vectorPtr_Type& u_star, const vectorPtr_Type& rhs_velocity );

	// Solve the current timestep, provided the BC
	void iterate( bcPtr_Type & bc, const Real& time );

	// Set coefficient associated to the time discretization scheme
	void setAlpha(const Real& alpha)
	{
		M_alpha = alpha;
	}

	void setParameters( );

	// Set coefficient associated to the time discretization scheme
	void setTimeStep(const Real& dt)
	{
		M_timeStep = dt;
	}

	// Get the velocity FE space
	const boost::shared_ptr<FESpace<mesh_Type, map_Type> >& uFESpace() const
	{
		return M_velocityFESpace;
	}

	// Get the velocity FE space
	const boost::shared_ptr<FESpace<mesh_Type, map_Type> >& pFESpace() const
	{
		return M_pressureFESpace;
	}

    void updateVelocity(vectorPtr_Type& velocity)
    {
        *velocity = *M_velocity;
    }
    
    void updatePressure(vectorPtr_Type& pressure)
    {
        *pressure = *M_pressure;
    }
    
    void applyBoundaryConditions ( bcPtr_Type & bc, const Real& time );
    
    matrixPtr_Type const& getF() const
    {
    	return M_F;
    }

    matrixPtr_Type const& getBtranspose() const
    {
    	return M_Btranspose;
    }

    matrixPtr_Type const& getB() const
    {
    	return M_B;
    }

    vectorPtr_Type const& getRhs() const
    {
    	return M_rhs;
    }

private:

	// build the graphs
	void buildGraphs();

	// update the bc handler
	void updateBCHandler( bcPtr_Type & bc );

    void setSolversOptions(const Teuchos::ParameterList& solversOptions);

	// communicator
	commPtr_Type M_comm;

	// getpot object
	dataFile_Type M_dataFile;

	// fluid data
	boost::shared_ptr<OseenData> M_fluidData;

	// FE spaces
	boost::shared_ptr<FESpace<mesh_Type, map_Type> > M_velocityFESpace;
	boost::shared_ptr<FESpace<mesh_Type, map_Type> > M_pressureFESpace;

	// ET FE Spaces
	boost::shared_ptr<ETFESpace_velocity > M_fespaceUETA;
	boost::shared_ptr<ETFESpace_pressure > M_fespacePETA;

	// stiff-strain check
	bool M_stiffStrain;

	// graphs for the matrices
	graphPtr_Type M_Mu_graph;
	graphPtr_Type M_Btranspose_graph;
	graphPtr_Type M_B_graph;
	graphPtr_Type M_C_graph;
	graphPtr_Type M_A_graph;
	graphPtr_Type M_F_graph;

	// matrices
	matrixPtr_Type M_Mu;
	matrixPtr_Type M_Btranspose;
	matrixPtr_Type M_B;
	matrixPtr_Type M_C;
	matrixPtr_Type M_A;
	matrixPtr_Type M_F;

	// vectors
	vectorPtr_Type M_uExtrapolated;
	vectorPtr_Type M_rhs;
    vectorPtr_Type M_velocity;
    vectorPtr_Type M_pressure;
    
	//! Displayer to print in parallel (only PID 0 will print)
	Displayer M_displayer;

	//! Reals
	Real M_alpha;
	Real M_timeStep;
	bool M_graphIsBuilt;
    
    // Navoer Stokes operator
	boost::shared_ptr<LifeV::Operators::NavierStokesOperator> M_oper;
    
    // Navoer Stokes operator
	boost::shared_ptr<LifeV::Operators::aSIMPLEOperator> M_prec;
    
    // Epetra Operator needed to solve the linear system
    boost::shared_ptr<Operators::InvertibleOperator> M_invOper;
    
    // Parameter list solver
    parameterListPtr_Type M_pListLinSolver;

}; // class NavierStokesSolver

} // namespace LifeV

#endif // NAVIERSTOKESSOLVERBLOCK_H
