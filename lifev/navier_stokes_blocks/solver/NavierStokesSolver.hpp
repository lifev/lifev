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
#include <lifev/navier_stokes_blocks/solver/NavierStokesOperator.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>

#include <lifev/navier_stokes_blocks/solver/NavierStokesPreconditionerOperator.hpp>
#include <lifev/navier_stokes_blocks/solver/aSIMPLEOperator.hpp>
#include <lifev/navier_stokes_blocks/solver/aPCDOperator.hpp>

// utilities
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/util/Displayer.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <lifev/core/algorithm/NonLinearRichardson.hpp>

#include <lifev/core/filter/GetPot.hpp>

#include <lifev/navier_stokes_blocks/solver/StabilizationSUPG.hpp>

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

	// Solve the steady Navier-Stokes equations, provided the BC
	void iterate_steady( bcPtr_Type & bc );

	// Solve the Navier-Stokes equations, provided the BC
	void iterate_nonlinear( bcPtr_Type & bc, const Real& time );

	void evalResidual(vector_Type& residual, const vector_Type& solution, const UInt iter_newton);

	void solveJac( vector_Type& increment, const vector_Type& residual, const Real linearRelTol );

	// Update the Jacobian matrix
	void updateJacobian( const vector_Type& u_k );

	// Apply the BCs on the Jacobian matrix
	void applyBoundaryConditionsJacobian ( bcPtr_Type & bc );

	void applyBoundaryConditionsSolution ( bcPtr_Type & bc, const Real& time );

	void applyGravityForce ( const Real& gravity, const Real& gravityDirection);

	void updateConvectiveTerm ( const vectorPtr_Type& velocity);

	//! Evaluates the fluid residual in FSI simulations
	/*!
	 * @param convective_velocity difference between fluid velocity and mesh velocity at previous Newton iteration
	 * @param velocity_km1 fluid velocity at previous Newton iteration
	 * @param pressure_km1 fluid pressure at previous Newton iteration
	 * @param pressure_km1 terms of bdf associated to the rhs
	 * @param residual residual that will be assembled
	 */
	void evaluateResidual( const vectorPtr_Type& convective_velocity,
						   const vectorPtr_Type& velocity_km1,
						   const vectorPtr_Type& pressure_km1,
						   const vectorPtr_Type& rhs_velocity,
						   vectorPtr_Type& residual);

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

    matrixPtr_Type const& getJacobian() const
    {
    	return M_Jacobian;
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

    void setBCpcd(const bcPtr_Type & bc)
    {
    	M_bcPCD = bc;
    }

    void updatePCD(const vectorPtr_Type& velocity);

    matrixPtr_Type const& Fp() const
    {
    	return M_Fp;
    }

    matrixPtr_Type const& Mp() const
    {
    	return M_Mp;
    }

    matrixPtr_Type const& Mu() const
    {
    	return M_Mu;
    }

    void applyBCMu(bcPtr_Type& pcdBC)
    {
    	bcManageMatrix( *M_Mu, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *pcdBC, M_velocityFESpace->feBd(), 1.0, 0.0);
    	M_Mu->globalAssemble();
    }

    void setBC ( const bcPtr_Type& bc )
    {
    	M_bc = bc;
    }

    Real density (  ) const
    {
    	return M_density;
    }

    Real viscosity (  ) const
    {
    	return M_viscosity;
    }

    void setRhsVelocity ( const vectorPtr_Type& vel_rhs)
    {
    	M_velocityRhs.reset( new vector_Type ( *vel_rhs ) );
    }

    matrixPtr_Type const& block00() const
    {
    	return M_block00;
    }

    matrixPtr_Type const& block01() const
    {
    	return M_block01;
    }

    matrixPtr_Type const& block10() const
    {
    	return M_block10;
    }

    matrixPtr_Type const& block11() const
    {
    	return M_block11;
    }

private:

	// build the graphs
	void buildGraphs();

	// update the bc handler
	void updateBCHandler( bcPtr_Type & bc );

    void setSolversOptions(const Teuchos::ParameterList& solversOptions);

    void buildPCDGraphs();

	// communicator
	commPtr_Type M_comm;

	// getpot object
	dataFile_Type M_dataFile;

	// FE spaces
	boost::shared_ptr<FESpace<mesh_Type, map_Type> > M_velocityFESpace;
	boost::shared_ptr<FESpace<mesh_Type, map_Type> > M_pressureFESpace;
	boost::shared_ptr<FESpace<mesh_Type, map_Type> > M_velocityFESpaceScalar;

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
	graphPtr_Type M_Jacobian_graph;

	graphPtr_Type M_Mp_graph;
	graphPtr_Type M_Fp_graph;

	// matrices
	matrixPtr_Type M_Mu;
	matrixPtr_Type M_Btranspose;
	matrixPtr_Type M_B;
	matrixPtr_Type M_C;
	matrixPtr_Type M_A;
	matrixPtr_Type M_F;
	matrixPtr_Type M_Jacobian;

	matrixPtr_Type M_Mp;
	matrixPtr_Type M_Fp;

	// Navier-Stokes block matrices
	matrixPtr_Type M_block00;
	matrixPtr_Type M_block01;
	matrixPtr_Type M_block10;
	matrixPtr_Type M_block11;

	// vectors
	vectorPtr_Type M_uExtrapolated;
	vectorPtr_Type M_rhs;
    vectorPtr_Type M_velocity;
    vectorPtr_Type M_pressure;
    vectorPtr_Type M_velocityRhs;

    boost::shared_ptr<map_Type> M_monolithicMap;
    vectorPtr_Type M_solution;

    vectorPtr_Type M_residual_u;
    vectorPtr_Type M_residual_p;

    vectorPtr_Type M_velocity_old_newton;
    vectorPtr_Type M_pressure_old_newton;

	//! Displayer to print in parallel (only PID 0 will print)
	Displayer M_displayer;

	//! Reals
	Real M_alpha;
	Real M_timeStep;

	//! Booleans
	bool M_graphIsBuilt;
	bool M_graphPCDisBuilt;

    // Navoer Stokes operator
	boost::shared_ptr<LifeV::Operators::NavierStokesOperator> M_oper;

    // Preconditioner
	boost::shared_ptr<LifeV::Operators::NavierStokesPreconditionerOperator> M_prec;

    // Epetra Operator needed to solve the linear system
    boost::shared_ptr<Operators::InvertibleOperator> M_invOper;

    // Parameter list solver
    parameterListPtr_Type M_pListLinSolver;

    // Check if the convective term is fully implicit
    bool M_fullyImplicit;

    // BC for the pcd preconditioner
    bcPtr_Type M_bcPCD;

    // Check if we want to compute the steady state
    bool M_steady;

    Real M_relativeTolerance, M_absoluteTolerance, M_etaMax;
    Int  M_nonLinearLineSearch;
    UInt M_maxiterNonlinear;
    std::ofstream M_out_res;

    // BC handler
    bcPtr_Type M_bc;

    Real M_density;
    Real M_viscosity;

}; // class NavierStokesSolver

} // namespace LifeV

#endif // NAVIERSTOKESSOLVERBLOCK_H
