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
#define NAVIERSTOKESSOLVERBLOCK_H 1

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
//#include <lifev/navier_stokes/solver/SolverManager.hpp>
//#include <lifev/navier_stokes/solver/solverManager_aSIMPLE.hpp>

// utilities
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/util/Displayer.hpp>

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
    
    //! Setter for the solvers options
    void setSolversOptions (const Teuchos::ParameterList& solversOptions);

private:

	// build the graphs
	void buildGraphs();

	// update the bc handler
	void updateBCHandler( bcPtr_Type & bc );

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

	//! Displayer to print in parallel (only PID 0 will print)
	Displayer M_displayer;

	//! Reals
	Real M_alpha;
	Real M_timeStep;

	// Solver Manager
	// boost::shared_ptr<SolverManager> M_solverManager;
    
    // Navoer Stokes operator
	boost::shared_ptr<LifeV::Operators::NavierStokesOperator> M_oper;
    
    // Epetra Operator needed to solve the linear system
    // invOpPtr_Type M_invOp;

}; // class NavierStokesSolver


NavierStokesSolver::NavierStokesSolver(const dataFile_Type dataFile, const commPtr_Type& communicator):
		M_comm(communicator),
		M_dataFile(dataFile),
		M_displayer(communicator),
        M_oper(new Operators::NavierStokesOperator)
{
}

NavierStokesSolver::~NavierStokesSolver()
{
}


void NavierStokesSolver::setup(const meshPtr_Type& mesh)
{
	M_fluidData.reset( new OseenData() );
	M_fluidData->setup( M_dataFile );

	std::string uOrder = M_dataFile("fluid/space_discretization/vel_order","P1");
	std::string pOrder = M_dataFile("fluid/space_discretization/pres_order","P1");

	Int geoDimensions = mesh_Type::S_geoDimensions;

	M_velocityFESpace.reset (new FESpace<mesh_Type, map_Type> (mesh, uOrder, geoDimensions, M_comm) );
	M_pressureFESpace.reset (new FESpace<mesh_Type, map_Type> (mesh, pOrder, geoDimensions, M_comm) );

	M_fespaceUETA.reset( new ETFESpace_velocity(M_velocityFESpace->mesh(), &(M_velocityFESpace->refFE()), M_comm));
	M_fespacePETA.reset( new ETFESpace_pressure(M_pressureFESpace->mesh(), &(M_pressureFESpace->refFE()), M_comm));

	M_uExtrapolated.reset( new vector_Type ( M_velocityFESpace->map(), Repeated ) );
	*M_uExtrapolated *= 0;

	M_stiffStrain = M_dataFile("fluid/space_discretization/stiff_strain", true);

	//std::string solverManagerType = M_dataFile("fluid/solver_manager_type"," ");
	//M_solverManager.reset(SolverFactory::instance().createObject(solverManagerType));
	//M_solverManager->setComm(M_comm);
}

void NavierStokesSolver::buildGraphs()
{
	M_displayer.leaderPrint ( " F - Pre-building the graphs... ");
	LifeChrono chrono;
	chrono.start();

	{
		using namespace ExpressionAssembly;

		// Graph velocity mass -> block (0,0)
		M_Mu_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_fespaceUETA->mesh() ),
					 quadRuleTetra4pt,
					 M_fespaceUETA,
					 M_fespaceUETA,
					 M_fluidData->density() * dot ( phi_i, phi_j )
		) >> M_Mu_graph;
		M_Mu_graph->GlobalAssemble();
        M_Mu_graph->OptimizeStorage();

		// Graph block (0,1) of NS
		M_Btranspose_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_fespaceUETA->mesh() ),
					 quadRuleTetra4pt,
					 M_fespaceUETA,
					 M_fespacePETA,
					 value(-1.0) * phi_j * div(phi_i)
		) >> M_Btranspose_graph;
		M_Btranspose_graph->GlobalAssemble( *(M_pressureFESpace->map().map (Unique)), *(M_velocityFESpace->map().map (Unique)) );
        M_Btranspose_graph->OptimizeStorage();
        
		// Graph block (1,0) of NS
		M_B_graph.reset (new Epetra_FECrsGraph (Copy, *(M_pressureFESpace->map().map (Unique)), 0) );
		buildGraph ( elements (M_fespaceUETA->mesh() ),
					 quadRuleTetra4pt,
					 M_fespacePETA,
					 M_fespaceUETA,
					 phi_i * div(phi_j)
		) >> M_B_graph;
		M_B_graph->GlobalAssemble( *(M_velocityFESpace->map().map (Unique)), *(M_pressureFESpace->map().map (Unique)) );
        M_B_graph->OptimizeStorage();
        
		// Graph convective term, block (0,0)
		M_C_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_fespaceUETA->mesh() ),
					 quadRuleTetra4pt,
					 M_fespaceUETA,
					 M_fespaceUETA,
					 dot( M_fluidData->density()*value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i)
		) >> M_C_graph;
		M_C_graph->GlobalAssemble();
        M_C_graph->OptimizeStorage();

		// Graph stiffness, block (0,0)
		M_A_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
		if (M_stiffStrain)
		{
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 value( 0.5 * M_fluidData->viscosity() ) * dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) )
			) >> M_A_graph;
		}
		else
		{
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 M_fluidData->viscosity() * dot( grad(phi_i) , grad(phi_j) + transpose(grad(phi_j)) )
			) >> M_A_graph;
		}
		M_A_graph->GlobalAssemble();
        M_A_graph->OptimizeStorage();
        
		// Graph of entire block (0,0)
		M_F_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
		if (M_stiffStrain)
		{
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 dot ( phi_i, phi_j ) + // mass
						 dot( M_fluidData->density()*value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i) + // convective term
						 value( 0.5 * M_fluidData->viscosity() ) * dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) ) // stiffness
			) >> M_F_graph;
		}
		else
		{
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 dot ( phi_i, phi_j ) + // mass
						 dot( M_fluidData->density()*value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i) + // convective term
						 M_fluidData->viscosity() * dot( grad(phi_i) , grad(phi_j) + transpose(grad(phi_j)) ) // stiffness
			) >> M_F_graph;
		}
		M_F_graph->GlobalAssemble();
        M_F_graph->OptimizeStorage();
	}

	chrono.stop();
	M_displayer.leaderPrintMax ( "   done in ", chrono.diff() ) ;
}

void NavierStokesSolver::buildSystem()
{
	buildGraphs();

	M_displayer.leaderPrint ( " F - Assembling constant terms... ");
	LifeChrono chrono;
	chrono.start();

	{
		using namespace ExpressionAssembly;

		// Graph velocity mass -> block (0,0)
		M_Mu.reset (new matrix_Type ( M_velocityFESpace->map(), *M_Mu_graph ) );
		*M_Mu *= 0;
		integrate ( elements (M_fespaceUETA->mesh() ),
					M_velocityFESpace->qr(),
					M_fespaceUETA,
					M_fespaceUETA,
					M_fluidData->density() * dot ( phi_i, phi_j )
		) >> M_Mu;
		M_Mu->globalAssemble();

		M_Btranspose.reset (new matrix_Type ( M_velocityFESpace->map(), *M_Btranspose_graph ) );
		*M_Btranspose *= 0;
		integrate( elements(M_fespaceUETA->mesh()),
				   M_velocityFESpace->qr(),
				   M_fespaceUETA,
				   M_fespacePETA,
				   value(-1.0) * phi_j * div(phi_i)
		) >> M_Btranspose;
		M_Btranspose->globalAssemble( M_pressureFESpace->mapPtr(), M_velocityFESpace->mapPtr() );

		M_B.reset (new matrix_Type ( M_pressureFESpace->map(), *M_B_graph ) );
		*M_B *= 0;
		integrate( elements(M_fespaceUETA->mesh()),
				   M_pressureFESpace->qr(),
				   M_fespacePETA,
				   M_fespaceUETA,
				   phi_i * div(phi_j)
		) >> M_B;
		M_B->globalAssemble( M_velocityFESpace->mapPtr(), M_pressureFESpace->mapPtr());

		M_A.reset (new matrix_Type ( M_velocityFESpace->map(), *M_A_graph ) );
		*M_A *= 0;
		if ( M_stiffStrain )
		{
			integrate( elements(M_fespaceUETA->mesh()),
					   M_velocityFESpace->qr(),
					   M_fespaceUETA,
					   M_fespaceUETA,
					   value( 0.5 * M_fluidData->viscosity() ) * dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) )
			) >> M_A;
		}
		else
		{
			integrate( elements(M_fespaceUETA->mesh()),
					   M_velocityFESpace->qr(),
					   M_fespaceUETA,
				 	   M_fespaceUETA,
					   M_fluidData->viscosity() * dot( grad(phi_i) , grad(phi_j) + transpose(grad(phi_j)) )
			) >> M_A;
		}
		M_A->globalAssemble();

		M_C.reset (new matrix_Type ( M_velocityFESpace->map(), *M_C_graph ) );
		*M_C *= 0;
		M_C->globalAssemble();

		M_F.reset (new matrix_Type ( M_velocityFESpace->map(), *M_F_graph ) );
		*M_F *= 0;
		M_F->globalAssemble();
	}

	chrono.stop();
	M_displayer.leaderPrintMax ( " done in ", chrono.diff() ) ;
}

void NavierStokesSolver::updateSystem( const vectorPtr_Type& u_star, const vectorPtr_Type& rhs_velocity )
{
	M_uExtrapolated.reset( new vector_Type ( *u_star, Repeated ) );

	// Update convective term
	*M_C *= 0;
	{
		using namespace ExpressionAssembly;
		integrate( elements(M_fespaceUETA->mesh()),
				   M_velocityFESpace->qr(),
				   M_fespaceUETA,
				   M_fespaceUETA,
				   dot( M_fluidData->density() * value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i) // semi-implicit treatment of the convective term
		)
		>> M_C;
	}
	M_C->globalAssemble();

	// Get the matrix corresponding to the block (0,0)
	*M_F *= 0;
	*M_F += *M_Mu;
	*M_F *= M_alpha/M_timeStep;
	*M_F += *M_A;
	*M_F += *M_C;

	// Get the right hand side with inertia contribution
	M_rhs.reset( new vector_Type ( M_velocityFESpace->map(), Repeated ) );
	*M_rhs = *M_Mu* (*rhs_velocity);
}

void NavierStokesSolver::iterate( bcPtr_Type & bc, const Real& time )
{
	updateBCHandler(bc);
	bcManage ( *M_F, *M_rhs, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 1.0, time );

	for (BCHandler::bcBaseIterator_Type it = bc->begin(); it != bc->end(); ++it)
	{
		if ( it->type() != Essential)
			continue;

		bcManageMatrix( *M_Btranspose, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 0.0, 0.0);
	}
    
    //(2) Set up the OseenOperator
    M_displayer.leaderPrint( "\tset up the block operator...");
    LifeChrono chrono;
    chrono.start();
    
    Operators::NavierStokesOperator::operatorPtrContainer_Type operData(2,2);
    operData(0,0) = M_F->matrixPtr();
    operData(0,1) = M_Btranspose->matrixPtr();
    operData(1,0) = M_B->matrixPtr();
    
    M_oper->setUp(operData, M_displayer.comm());
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );
    
}

void NavierStokesSolver::updateBCHandler( bcPtr_Type & bc )
{
	bc->bcUpdate ( *M_velocityFESpace->mesh(), M_velocityFESpace->feBd(), M_velocityFESpace->dof() );
}
    
//! Setter for the solvers options
void NavierStokesSolver::setSolversOptions (const Teuchos::ParameterList& solversOptions)
{
    /*
    boost::shared_ptr<Teuchos::ParameterList> schurComplementList;
    schurComplementList.reset(new Teuchos::ParameterList(solversOptions.sublist("ApproximatedSchurOperator")) );
    M_solverManager->setSchurOptions(schurComplementList);
    
    boost::shared_ptr<Teuchos::ParameterList> momentumList;
    momentumList.reset(new Teuchos::ParameterList(solversOptions.sublist("MomentumOperator")) );
    M_solverManager->setMomentumOptions(momentumList);
    
    boost::shared_ptr<Teuchos::ParameterList> monolithicList;
    monolithicList.reset(new Teuchos::ParameterList(solversOptions.sublist("MonolithicOperator")) );
    M_solverManager->setLinSolverParameter(monolithicList);
    */
}

} // namespace LifeV

#endif // NAVIERSTOKESSOLVERBLOCK_H
