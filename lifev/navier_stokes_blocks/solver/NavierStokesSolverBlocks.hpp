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
   @brief Navier Stokes solver.
   @author Davide Forti <davide.forti@epfl.ch>
   @maintainer Davide Forti <davide.forti@epfl.ch>
   @date 03-02-2015

   This class implements a Navier-Stokes solver. This solver allows to solve this kind of problems:
   - Steady NS equations (using P2-P1 or P1Bubble-P1 finite elements)
   - Unsteady NS equations using a fully-implicit time discretization (SUPG-VMS stabilization available)
   - Unsteady NS equations using a semi-implicit time discretization (SUPG-VMS and VMS-LES stabilizations available)

   In the testsuite folder several examples which show how to properly use the solver are present.

   If you are using this solver to generate results for publication, please <b>cite</b>:
   - D. Forti, L. Dede'. <i> Semi-implicit BDF time discretization of the Navier–Stokes equations with VMS-LES modeling in a High Performance Computing framework</i>.
   	 Comput. Fluids. 197(1):168-182.
 */

#ifndef NAVIERSTOKESSOLVERBLOCKS_H
#define NAVIERSTOKESSOLVERBLOCKS_H

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

// utilities
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/util/Displayer.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <lifev/core/algorithm/NonLinearRichardson.hpp>

#include <lifev/core/filter/GetPot.hpp>

#include <lifev/navier_stokes_blocks/solver/Stabilization.hpp>
#include <lifev/navier_stokes_blocks/solver/StabilizationSUPG.hpp>
#include <lifev/navier_stokes_blocks/solver/StabilizationSUPG_semi_implicit.hpp>
#include <lifev/navier_stokes_blocks/solver/StabilizationSUPGALE.hpp>
#include <lifev/navier_stokes_blocks/solver/StabilizationSUPG_semi_implicit_ale.hpp>
#include <lifev/navier_stokes_blocks/solver/StabilizationVMSLES_semi_implicit.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>

#include <lifev/core/fem/PostProcessingBoundary.hpp>
#include <lifev/navier_stokes_blocks/solver/FastAssemblerNS.hpp>

namespace LifeV
{

class NavierStokesSolverBlocks
{

public:

	//@name Public Types
	//@{

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

    typedef LifeV::Preconditioner                  basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>       basePrecPtr_Type;
    typedef LifeV::PreconditionerIfpack            prec_Type;
    typedef boost::shared_ptr<prec_Type>           precPtr_Type;
    typedef Teuchos::RCP< Teuchos::ParameterList > parameterListRCP_Type;

    //@}

    //! @name Constructor and Destructor
    //@{

	// Constructor
    NavierStokesSolverBlocks(const dataFile_Type dataFile, const commPtr_Type& communicator);

	// Destructor
	~NavierStokesSolverBlocks();

	//@}

	//! @name Methods
	//@{

	//! Setup for the solver
	/*!
	 * @param mesh mesh used
	 * @param id_domain used just in case of multiple fluids
	 */
	void setup(const meshPtr_Type& mesh, const int& id_domain = 36);

	//! Assemble constant terms
	void buildSystem();

	//! Update the convective term and the right hand side
	/*!
	 * @param u_star extrapolate velocity
	 * @param rhs_velocity right hand side, velocity block
	 */
	void updateSystem( const vectorPtr_Type& u_star, const vectorPtr_Type& rhs_velocity );

	//! Solve the current timestep, provided the BC
	/*!
	 * @param bc boundary conditions
	 * @param time current time
	 */
	void iterate( bcPtr_Type & bc, const Real& time );

	//! Solve the current timestep, provided the BC and a vector of velocities (used only for for aorta example)
	/*!
	  * @param bc boundary conditions
	  * @param time current time
	  * @param velocities vector of velocity (used for imposing velocity on noncircular inflows)
	 */
	void iterate( bcPtr_Type & bc, const Real& time, const vectorPtr_Type& velocities );

	//! Solve the steady Navier-Stokes equations
	void iterate_steady( );

	//! Solve the Navier-Stokes equations at a certain time
	/*!
	 * @param time current time
	 */
	void iterate_nonlinear( const Real& time );

	//! Evaluate the residual of the NS problem
	/*!
	 * @param residual residual vector
	 * @param solution solution vector
	 * @param iter_newton current newton iteration
	 */
	void evalResidual(vector_Type& residual, const vector_Type& solution, const UInt iter_newton);

	//! Evaluate the residual of the NS problem
	/*!
	 * @param increment increment vector
	 * @param residual residual vector
	 * @param linearRelTol tolerance of the linar solver
	 */
	void solveJac( vector_Type& increment, const vector_Type& residual, const Real linearRelTol );

	//! Update the Jacobian matrix, only the term associated with the linearization of the convective term
	/*!
	 * @param u_k velocity previous newton step
	 */
	void updateJacobian( const vector_Type& u_k );

	//! Apply the BCs semi-implicit case
	/*!
	 * @param bc boundary conditions of the problem
	 * @param time current time
	 */
    void applyBoundaryConditions ( bcPtr_Type & bc, const Real& time );

    //! Apply the BCs semi-implicit case (example aorta)
    /*!
     * @param bc boundary conditions of the problem
     * @param time current time
     * @param velocities vector of velocities to be imposed
     */
    void applyBoundaryConditions ( bcPtr_Type & bc, const Real& time, const vectorPtr_Type& velocities );

	//! Apply the BCs on the Jacobian matrix
	/*!
	 * @param bc boundary conditions of the problem
	 */
	void applyBoundaryConditionsJacobian ( bcPtr_Type & bc );

	//! Apply the BCs to the solution
	/*!
	 * @param time current time
	 */
	void applyBoundaryConditionsSolution ( const Real& time );

	//! Update the convective term
	/*!
	 *	 @param velocity velocity vector
	 */
	void updateConvectiveTerm ( const vectorPtr_Type& velocity);

	//! Compute Aerodynamic Forces
	/*!
	 *	 @param bcHandlerDrag bc to be used for drag computation
	 *	 @param bcHandlerLift bc to be used for lift computation
	 */
    VectorSmall<2> computeForces( BCHandler& bcHandlerDrag, BCHandler& bcHandlerLift);

    //! Solve a time step
    void solveTimeStep();

    //! Compute Aerodynamic Forces - nonlinear case
    /*!
   	 *	 @param force vector containing the forces
   	 *	 @param solution vector with the solution
   	 */
    void computeForcesNonLinear(vectorPtr_Type& force, const vectorPtr_Type& solution);

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

	//! Evaluates the fluid residual in FSI simulations
	/*!
	 * @param convective_velocity difference between fluid velocity and mesh velocity at previous Newton iteration
	 * @param velocity_km1 fluid velocity at previous Newton iteration
	 * @param pressure_km1 fluid pressure at previous Newton iteration
	 * @param pressure_km1 terms of bdf associated to the rhs
	 * @param residualVelocity residual associated to the fluid momentum equation
	 * @param residualPressure residual associated to the fluid continuity equation
	 */
	void evaluateResidual( const vectorPtr_Type& convective_velocity,
						   const vectorPtr_Type& velocity_km1,
						   const vectorPtr_Type& pressure_km1,
						   const vectorPtr_Type& rhs_velocity,
						   vectorPtr_Type& residualVelocity,
						   vectorPtr_Type& residualPressure);

	//! Update stabilization terms
	/*!
	 * @param convective_velocity_previous_newton_step difference between fluid velocity and mesh velocity at previous Newton iteration
	 * @param velocity_previous_newton_step fluid velocity at previous Newton iteration
	 * @param pressure_previous_newton_step fluid pressure at previous Newton iteration
	 * @param velocity_rhs block of the velocity of the rhs
	 */
	void updateStabilization( const vector_Type& convective_velocity_previous_newton_step,
							  const vector_Type& velocity_previous_newton_step,
							  const vector_Type& pressure_previous_newton_step,
							  const vector_Type& velocity_rhs );

	//! Computation of forces
	/*!
	 * @param velocity velocity vector
	 * @param pressure pressure vector
	 */
	void integrateForces ( const vectorPtr_Type & velocity, const vectorPtr_Type & pressure);

	//! Additional method used for pre-processing on non-circular boundaries (used only in example aorta)
	/*!
	 * @param nx x component normal vector to an outflow face
	 * @param ny y component normal vector to an outflow face
	 * @param nz z component normal vector to an outflow face
	 * @param bc boundary conditions
	 * @param Q_hat reference flowrate to be imposed
	 * @param Phi_h solution laplacian
	 * @param flag flag outflow
	 * @param Phi_h_flag solution laplacian on the outflow flag
	 * @param V_hat_x reference velocity, x component
	 * @param V_hat_y reference velocity, y component
	 * @param V_hat_z reference velocity, z component
	 */
	void preprocessBoundary(const Real& nx, const Real& ny, const Real& nz, BCHandler& bc, Real& Q_hat, const vectorPtr_Type& Phi_h, const UInt flag,
			vectorPtr_Type& Phi_h_flag, vectorPtr_Type& V_hat_x, vectorPtr_Type& V_hat_y, vectorPtr_Type& V_hat_z);

	//! Additional method used for pre-processing on non-circular boundaries (used only in example aorta)
	/*!
	 * @param flag flag of the outflow face
	 * @param bc_laplacian boundary conditions laplacian
	 * @param laplacianSolution solution laplacian
	 */
	void solveLaplacian( const UInt& flag, bcPtr_Type& bc_laplacian, vectorPtr_Type& laplacianSolution );

	//! Additional setup for postprocessing on boundaries
	void setupPostProc( );

	//! Compute flux through a boundary face of the domain
	/*!
	 * @param flag flag of the face
	 * @param velocity velocity vector
	 */
	Real flux ( const markerID_Type& flag, const vector_Type& velocity );

	//! Compute area of a boundary face of the domain
	/*!
	 * @param flag flag of the outflow face
	 */
	Real area ( const markerID_Type& flag );

	//! Compute mean pressure at a boundary face of the domain
	/*!
	 * @param flag flag of the outflow face
	 * @param pressure pressure vector
	 */
	Real pres ( const markerID_Type& flag, const vector_Type& pressure );

	//! Compute center of a boundary face of the domain
	/*!
	 * @param flag flag of the face
	 */
	Vector geometricCenter ( const markerID_Type& flag );

	//! Compute normal of a boundary face of the domain
	/*!
	 * @param flag flag of the face
	 */
	Vector normal ( const markerID_Type& flag );

	//@}

	//! @name Set methods
	//@{

	//! Set coefficient associated to the time discretization scheme
	/*!
	 *	@param alpha coefficient BDF scheme in front of u^{n+1}
	 */
	void setAlpha(const Real& alpha)
	{
		M_alpha = alpha;
		if ( M_useStabilization )
		{
			M_stabilization->setAlpha( M_alpha );
		}
	}

	//! Set the rhs vector, velocity component
	/*!
	 *	@param vel_rhs rhs vector, velocity component
	 */
	void setRhsVelocity ( const vectorPtr_Type& vel_rhs)
	{
		M_velocityRhs.reset( new vector_Type ( *vel_rhs ) );
	}

	//! Set the extrapolated pressure vector (used by semi-implicit VMS-LES stabilization)
	/*!
	 * @param pressure_extrapolated pressure extrapolated vector
	 */
	void setExtrapolatedPressure( const vectorPtr_Type& pressure_extrapolated ) { M_pressure_extrapolated = pressure_extrapolated; }

	//! Set if one needs to export the fine scale component (when using LES models)
	/*!
	 *	@param exporter exporter
	 *	@param numElementsTotal number of total elements (tetrahedral elements)
	 */
	void setExportFineScaleVelocity(ExporterHDF5<mesh_Type>& exporter, const int& numElementsTotal);

	//! Set parameters of the solver
	void setParameters( );

	//! Set time step
	/*!
	 *	@param dt time step size
	 */
	void setTimeStep(const Real& dt)
	{
		M_timeStep = dt;
	}

	//! Set time step
	/*!
	 *	@param bc boundary conditions of the problem
	 */
	void setBoundaryConditions ( const bcPtr_Type& bc );

	//@}

	//! @name Get methods
	//@{

	//! Get the velocity FE space
	const boost::shared_ptr<FESpace<mesh_Type, map_Type> >& uFESpace() const
	{
		return M_velocityFESpace;
	}

	//! Get the velocity FE space
	const boost::shared_ptr<FESpace<mesh_Type, map_Type> >& pFESpace() const
	{
		return M_pressureFESpace;
	}

	//! Get the velocity FE space
	const boost::shared_ptr<FESpace<mesh_Type, map_Type> >& uFESpace_scalar() const
	{
		return M_velocityFESpaceScalar;
	}

	//! Get the velocity vector
    void updateVelocity(vectorPtr_Type& velocity)
    {
        *velocity = *M_velocity;
    }

    //! Get the pressure vector
    void updatePressure(vectorPtr_Type& pressure)
    {
        *pressure = *M_pressure;
    }

    //! Get the F block
    matrixPtr_Type const& getF() const
    {
    	return M_F;
    }

    //! Get the jacobian block
    matrixPtr_Type const& getJacobian() const
    {
    	return M_Jacobian;
    }

    //! Get the Btranspose block
    matrixPtr_Type const& getBtranspose() const
    {
    	return M_Btranspose;
    }

    //! Get the B block
    matrixPtr_Type const& getB() const
    {
    	return M_B;
    }

    //! Get the rhs
    vectorPtr_Type const& getRhs() const
    {
    	return M_rhs;
    }

    //! Get the rhs vector (pressure part)
    vectorPtr_Type const& getRhsPressure() const
    {
    	return M_rhs_pressure;
    }

    matrixPtr_Type const& Mu() const
    {
    	return M_Mu;
    }
    
    Real density (  ) const
    {
    	return M_density;
    }

    Real viscosity (  ) const
    {
    	return M_viscosity;
    }

    //! Get the (0,0) block
    matrixPtr_Type const& block00() const
    {
    	return M_block00;
    }

    //! Get the (0,1) block
    matrixPtr_Type const& block01() const
    {
    	return M_block01;
    }

    //! Get the (1,0) block
    matrixPtr_Type const& block10() const
    {
    	return M_block10;
    }

    //! Get the (1,1)
    matrixPtr_Type const& block11() const
    {
    	return M_block11;
    }

    //! Get if using a stabilization
    bool useStabilization() const
    {
    	return M_useStabilization;
    }

    void assembleInterfaceMass( matrixPtr_Type& mass_interface, const mapPtr_Type& interface_map,
    						    markerID_Type interfaceFlag, const vectorPtr_Type& numerationInterface, const UInt& offset );

    //! Get the displayer
    Displayer const& getDisplayer ( ) const { return M_displayer; }

    //! Get the (0,0) block without BCs
    matrixPtr_Type block00_noBC() const
    {
    	return M_block00_noBC;
    }

    //! Get the (0,1) block without BCs
    matrixPtr_Type block01_noBC() const
    {
    	return M_block01_noBC;
    }

    //! Get the rhs without bcs applied
    vectorPtr_Type rhs_noBC() const
    {
    	return M_rhs_noBC;
    }

    //! Get the forces
    vectorPtr_Type getForces() const
    {
    	return M_forces;
    }

    //@}
    
private:

	//! Build the graphs
	void buildGraphs();

	//! Update the bc handler
	void updateBCHandler( bcPtr_Type & bc );

	//! Set options preconditioner
    void setSolversOptions(const Teuchos::ParameterList& solversOptions);

    //! Apply Bc to the residual vector
    void applyBoundaryConditionsResidual ( vector_Type& r_u, const Real& time = 0.0);
    
	// communicator
	commPtr_Type M_comm;

	// getpot object
	dataFile_Type M_dataFile;

	// Order FE spaces
	std::string M_uOrder;
	std::string M_pOrder;

	// FE spaces
	boost::shared_ptr<FESpace<mesh_Type, map_Type> > M_velocityFESpace;
	boost::shared_ptr<FESpace<mesh_Type, map_Type> > M_pressureFESpace;
	boost::shared_ptr<FESpace<mesh_Type, map_Type> > M_velocityFESpaceScalar;

	// ET FE Spaces
	boost::shared_ptr<ETFESpace_velocity > M_fespaceUETA;
	boost::shared_ptr<ETFESpace_pressure > M_fespacePETA;
	boost::shared_ptr<ETFESpace_pressure > M_fespaceUETA_scalar;

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

	// Navier-Stokes block matrices without boundary conditions applied
	// used for estimation of drag and lift coefficients when needed
	matrixPtr_Type M_block00_noBC;
	matrixPtr_Type M_block01_noBC;
	matrixPtr_Type M_block10_noBC;
	matrixPtr_Type M_block11_noBC;
	vectorPtr_Type M_rhs_noBC;
	vectorPtr_Type M_forces;
	boost::shared_ptr<LifeV::Operators::NavierStokesOperator> M_operLoads;

	// vectors
	vectorPtr_Type M_uExtrapolated;
	vectorPtr_Type M_rhs;
    vectorPtr_Type M_rhs_pressure;
    vectorPtr_Type M_velocity;
    vectorPtr_Type M_pressure;
    vectorPtr_Type M_velocityRhs;
    vectorPtr_Type M_velocityExtrapolated;

    boost::shared_ptr<map_Type> M_monolithicMap;
    vectorPtr_Type M_solution;

    vectorPtr_Type M_residual_u;
    vectorPtr_Type M_residual_p;

    vectorPtr_Type M_velocity_old_newton;
    vectorPtr_Type M_pressure_old_newton;
    
    vectorPtr_Type M_pressure_extrapolated;

	//! Displayer to print in parallel (only PID 0 will print)
	Displayer M_displayer;

	//! Reals
	Real M_alpha;
	Real M_timeStep;
    Real M_time;
    
	//! Booleans
	bool M_useGraph;
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
    bcPtr_Type M_bc_sol;
    bcPtr_Type M_bc_res_essential;
    bcPtr_Type M_bc_res_natural;
    
    Real M_density;
    Real M_viscosity;

    boost::shared_ptr<Stabilization> M_stabilization;

    bool M_useStabilization;
    
    std::string M_stabilizationType;

    bool M_nonconforming;

    // Members for imposition of weak boundary conditions
    bool M_imposeWeakBC;
    UInt M_flagWeakBC;
    Real M_meshSizeWeakBC;

    // Members for computation of aerodynamic loads
	bool M_computeAerodynamicLoads;
	std::string M_methodAerodynamicLoads;
	UInt M_flagBody;

	bool M_penalizeReverseFlow;
	UInt M_flagPenalizeReverseFlow;

	bool M_solve_blocks;
    
    //! Postprocessing class
    boost::shared_ptr<PostProcessingBoundary<mesh_Type> > M_postProcessing;

    bool M_useFastAssembly;
    boost::shared_ptr<FastAssemblerNS> M_fastAssembler;
    Real M_orderBDF;
    UInt M_orderVel;

}; // class NavierStokesSolverBlocks

} // namespace LifeV

#endif // NAVIERSTOKESSOLVERBLOCK_H
