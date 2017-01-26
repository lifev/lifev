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
@brief FSIHandler - File handling the solution of the FSI problem

@author Davide Forti <davide.forti@epfl.ch>
@date 28-10-2014

@maintainer Antonello Gerbi <antonello.gerbi@epfl.ch>
 
 This class implements a monolithic fully-implicit FSI solver. This solver features:
 - ALE treatment of the moving fluid domain
 - Linear/Nonlinear solid models can be used
 - Conforming/Nonconforming discretizations can be used
 - For the solid time discretizations one may use BDF of Newmark schemes
 - FaCSI preconitioner
 - In the fluid problem, if one chooses P1-P1 FE, a SUPG stabilization stemming from VMS approach can be used.
 - Mesh motion algorithms based on harmonic extensions or linear elasticity.
 - Multilayered structural model: a thin layer is considered underneath the thick 3D one. The thin layer is modeled 
 by a membrane model.
 
 In the testsuite folder, several examples which show how to properly use the solver are present.
 
 If you are using this solver to generate results for publication, please <b>cite</b>:
 - S. Deparis, D. Forti, G. Grandperrin, A. Quarteroni. <i> FaCSI: a block parallel preconditioner for fluid
 structure interaction in hemodynamics</i>. JCP (327), 700-718, 2016.
 - D. Forti, L. Dede'. <i> Semi-implicit BDF time discretization of the Navier–Stokes equations with VMS-LES modeling in a High Performance 
 Computing framework</i>. Comput. Fluids. 197(1):168-182.
 - S. Deparis, D. Forti, A. Quarteroni. <i> A fluid-structure interaction algorithm using radial basis function interpolation between non
 conforming interfaces</i>. In Y. Bazilevs, K. Takizawa (eds.) Advances in Computational Fluid-Structure, Modeling and Simulation in Science,
 Engineering and Technology. Birkhäuser Basel, 2016.
 - S. Deparis, D. Forti, P. Gervasio, and A. Quarteroni. <i>INTERNODES: an accurate interpolation-based method for coupling the Galerkin 
 solutions of PDEs on subdomains featuring non-conforming interfaces</i>. Comput. Fluids. 2016
 - D. Forti, M. Bukac, A. Quaini, S. Canic, S. Deparis. A monolithic approach to solving fluid-structure interaction between multilayered
 structures and incompressible, viscous fluids. Accepted in JSC, 2017.
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
#include <lifev/navier_stokes_blocks/solver/NavierStokesSolverBlocks.hpp>
#include <lifev/fsi_blocks/solver/LinearElasticity.hpp>

#include <lifev/fsi_blocks/solver/ALESolver.hpp>

// Expression template FE space
#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/core/fem/TimeAndExtrapolationHandler.hpp>
#include <lifev/core/fem/Newmark.hpp>
#include <lifev/core/fem/BDFSecondOrderDerivative.hpp>

#include <lifev/core/fem/DOFInterface3Dto3D.hpp>

#include <lifev/fsi_blocks/solver/FSIcouplingCE.hpp>
#include <lifev/fsi_blocks/solver/FSIApplyOperator.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <lifev/fsi_blocks/solver/DirichletNeumannPreconditioner.hpp>

#include <lifev/core/filter/PartitionIO.hpp>
#include <lifev/fsi_blocks/filter/DOFInterfaceIO.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <lifev/core/interpolation/Interpolation.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>

#include <lifev/fsi_blocks/solver/FSIApplyOperatorNonConforming.hpp>

#include <lifev/fsi_blocks/solver/NeoHookean.hpp>

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

    typedef BCHandler bc_Type;
    typedef boost::shared_ptr<BCHandler> bcPtr_Type;

    typedef VectorEpetra vector_Type;
	typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

	typedef MatrixEpetra<Real> matrix_Type;
	typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

	typedef Teuchos::ParameterList parameterList_Type;
	typedef boost::shared_ptr<parameterList_Type> parameterListPtr_Type;

	typedef Interpolation interpolation_Type;
	typedef boost::shared_ptr<interpolation_Type> interpolationPtr_Type;

	typedef LifeV::Preconditioner            basePrec_Type;
	typedef boost::shared_ptr<basePrec_Type> basePrecPtr_Type;

	typedef LifeV::PreconditionerIfpack  prec_Type;
	typedef boost::shared_ptr<prec_Type> precPtr_Type;

	typedef Teuchos::RCP< Teuchos::ParameterList > parameterListRCP_Type;

    //! Constructor
    FSIHandler(const commPtr_Type& communicator);

    //! Destructor
    ~FSIHandler();

    //! Set the datafile.
    /*!
     @param dataFile
     */
    void setDatafile(const GetPot& dataFile);

    //! Read the fluid and solid meshes.
    void readMeshes( );

    //! Partitioning the fluid and solid meshes.
    void partitionMeshes ( );

    //! Read fluid and solid meshes that have been already partitioned offline
    void readPartitionedMeshes( );

    //! Setup the solver
    void setup ( );

    //! Set the boundary conditions
    /*!
     @param fluidBC boundary conditions fluid problem
     @param structureBC boundary conditions solid problem
     @param aleBC boundary conditions ALE problem
     */
    void setBoundaryConditions ( const bcPtr_Type& fluidBC, const bcPtr_Type& structureBC, const bcPtr_Type& aleBC );

    //! Set the boundary conditions
    /*!
     @param interfaceFluidBC set dirichlet BC at the fluid interface for prec.
     */
    void setFluidInterfaceBoundaryConditions ( const bcPtr_Type& interfaceFluidBC ) { M_interfaceFluidBC = interfaceFluidBC; };

    //! Update all the bc handlers
    void updateBoundaryConditions( );

    //! Update the time advancing schemes
    void initializeTimeAdvance ( );

    //! Build the interface map
    void buildInterfaceMaps ( );

    //! Assemble the coupling blocks
    void assembleCoupling ( );

    //! Solves the time loop of the FSI problem
    void solveFSIproblem();

    //! Form the residual of the FSI problem
    /*!
     @param residual residual of the FSI problem
     @param solution current solution of FSI at certain Newton iteration
     @param iter_newton iteration Newton
     */
    void evalResidual(vector_Type& residual, const vector_Type& solution, const UInt iter_newton);

    //! Solves J_{FSI} \delta = - R
    /*!
     @param increment increment computed by Newton
     @param residual Rhs linear problem, residual FSI
     @param linearRelTol tolerance linear solver
     */
    void solveJac( vector_Type& increment, const vector_Type& residual, const Real linearRelTol );

    //! Set the parameter list of the problem
    void setParameterLists( );

    //! Set gravity, if considered
    /*!
     @param gravity value of the gravity
     @param gravity_direction 0 if x, 1 for y, 2 for z
    */
    void setGravity ( const Real& gravity, const Real& gravity_direction);

    //! Get the current time
    Real getTime () { return M_time; };

    //! Get initial time
    Real getStartTime () { return M_t_zero; };

    //! Set the current time
    void setTime (const Real& time ) { M_time = time; };

    //! Get the time step size
    Real getTimeStep () { return M_dt; };

    //! Get the end time
    Real getEndTime () { return M_t_end; };

    //! Get the fluid solver
    boost::shared_ptr<NavierStokesSolverBlocks> getFluid() { return M_fluid; };

    //! Method to be used before solveTimeStep. If one does not use
    // the solveFSIproblem() method, the following methods have to called:
    // intializeTimeLoop, solveTimeStep, postprocessResults and finalizeExporters.
    void intializeTimeLoop ( );

    //! Solve one time step.
    void solveTimeStep( );

    //! Save the results for post-processing.
    /*!
     @param time_step_count index of the time step
     */
    void postprocessResults( const int& time_step_count );

    //! Close the exporters
    void finalizeExporters ( );
    
    //! Get the fluid velocity
    vectorPtr_Type getFluidVelocity() { return M_fluidVelocity; };

    //! Get the fluid pressure
    vectorPtr_Type getFluidPressure() { return M_fluidPressure; };

    //! Get the FSI solution
    vectorPtr_Type getFSIsolution() { return M_solution; };

//@}

private:

    //! Moves the fluid mesh
    /*!
     @param displacement displacement fluid mesh
     */
    void moveMesh ( const VectorEpetra& displacement );

    //! Build the monolithic map
    void buildMonolithicMap ( );

    //! Setup solid sub-problem
    void setupStructure ( );

    //! Create ALE FE spaces
    void createAleFESpace();

    //! Update the bc handler
    void updateBCHandler( bcPtr_Type & bc );

    //! Build std map with local to global indexes dofs interface
    void createInterfaceMaps ( std::map<ID, ID> const& locDofMap );

    //! Build epetra map dofs at the interface
    void constructInterfaceMap ( const std::map<ID, ID>& locDofMap, const UInt subdomainMaxId);

    //! Setup the exporters of fluid and solid
    void setupExporters( );

    //! Instantiate exporter
    /*!
     @param exporter exporter
     @param localMesh mesh
     @param outputFileName filename output file
     */
    void instantiateExporter( boost::shared_ptr< Exporter<mesh_Type > > & exporter,
			  	  	  	  	  const meshPtr_Type& localMesh,
			  	  	  	  	  const std::string& outputFileName );
    
    //! Updates term due to coupling of velocities
    void updateSystem ( );

    //! Operator apply for the Jacobian in the onforming case
    void initializeApplyOperatorJacobian ( );

    //! Get Jacobian matrix structure problem
    void getMatrixStructure ( );

    //! Get vector coupling velocities at interface
    void get_structure_coupling_velocities ( );

    //! Apply BCs solid problem
    void applyBCstructure ( );

    //! Apply BCs solution FSI problem
    void applyBCsolution(vectorPtr_Type& M_solution);

    //! Apply BCs residual FSI problem
    /*!
     @param r_u residual momentum fluid
     @param r_ds residual momentum solid
     @param r_df residual ALE
     */
    void applyBCresidual(VectorEpetra& r_u, VectorEpetra& r_ds, VectorEpetra& r_df);

    //! Set options linear solver
    void setSolversOptions(const Teuchos::ParameterList& solversOptions);

    //! Update coupling velocity vector in the conforming case
    void updateRhsCouplingVelocities ( );

    //! Update coupling velocity vector in the nonconforming case
    void updateRhsCouplingVelocities_nonconforming ( );

    //! Extract interface dofs from vector defined on the whole domain
    void structureToInterface (vector_Type& VectorOnGamma, const vector_Type& VectorOnStructure);

    //! Initialize extrapolation for initial guess Newton
    void initializeExtrapolation( );

    //! Assemble interface mass structure at the interface
    void assembleStructureInterfaceMass ( );

    //! Apply inverse of interface mass matrix fluid
    /*!
     @param lambda residual weak form
     @param strongLambda residual strong form
    */
    void applyInverseFluidMassOnGamma ( const vectorPtr_Type& lambda, vectorPtr_Type& strongLambda );

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
    FESpacePtr_Type M_displacementFESpaceSerial;
    FESpacePtr_Type M_displacementFESpaceScalar;
    FESpacePtr_Type M_aleFESpace;

	solidETFESpacePtr_Type M_displacementETFESpace;

    // navier-stokes solver
    boost::shared_ptr<NavierStokesSolverBlocks> M_fluid;
    boost::shared_ptr<LinearElasticity> M_structure;
    boost::shared_ptr<NeoHookean> M_structureNeoHookean;
    boost::shared_ptr<ALESolver> M_ale;

    // time advance for the structure
    boost::shared_ptr<Newmark> M_structureTimeAdvance;
    boost::shared_ptr<TimeAndExtrapolationHandler> M_aleTimeAdvance;
    boost::shared_ptr<TimeAndExtrapolationHandler> M_fluidTimeAdvance;

    // boundary conditions
    bcPtr_Type M_fluidBC;
    bcPtr_Type M_fluidBC_residual_essential;
    bcPtr_Type M_fluidBC_residual_natural;
    bcPtr_Type M_structureBC;
    bcPtr_Type M_structureBC_residual_natural;
    bcPtr_Type M_structureBC_residual_essential;
    bcPtr_Type M_aleBC;
    bcPtr_Type M_aleBC_residual_natural;
    bcPtr_Type M_aleBC_residual_essential;
    bcPtr_Type M_interfaceFluidBC;
    

	//! Displayer to print in parallel (only PID 0 will print)
	Displayer M_displayer;

	//! Variables for the time advancing
	Real M_dt, M_t_zero, M_t_end, M_time;
	UInt M_orderBDF;

	//! Variables for the time advancing
	Real M_relativeTolerance, M_absoluteTolerance, M_etaMax;
	Int  M_nonLinearLineSearch;
	UInt M_maxiterNonlinear;
	std::ofstream M_out_res;

	boost::shared_ptr<DOFInterface3Dto3D> M_dofStructureToFluid;
	boost::shared_ptr<map_Type> M_structureInterfaceMap;
	boost::shared_ptr<map_Type> M_fluidInterfaceMap;
	boost::shared_ptr<map_Type> M_lagrangeMap;
	vectorPtr_Type M_numerationInterface;

	boost::shared_ptr<FSIcouplingCE> M_coupling;

	boost::shared_ptr<Exporter<mesh_Type > > M_exporterFluid;
	boost::shared_ptr<Exporter<mesh_Type > > M_exporterStructure;

	// Vectors for the exporters
	vectorPtr_Type M_fluidVelocity;
	vectorPtr_Type M_fluidPressure;
	vectorPtr_Type M_fluidDisplacement;
    vectorPtr_Type M_Lagrange;
    vectorPtr_Type M_LagrangeRestart;
	vectorPtr_Type M_structureDisplacement;
	vectorPtr_Type M_structureVelocity;
	vectorPtr_Type M_structureAcceleration;

	// Monolithic map
	boost::shared_ptr<map_Type> M_monolithicMap;

	vectorPtr_Type M_rhsFluid;
	vectorPtr_Type M_rhsStructure;
	vectorPtr_Type M_rhsCouplingVelocities;
	vectorPtr_Type M_solution;

	vectorPtr_Type M_beta;
	vectorPtr_Type M_rhs_velocity;

	matrixPtr_Type M_matrixStructure;
	matrixPtr_Type M_matrixStructure_noBc;

	// Operator to apply Jacobian of the FSI
	boost::shared_ptr<LifeV::Operators::FSIApplyOperator> M_applyOperatorJacobian;

	// Operator to compute the residual of the FSI
	boost::shared_ptr<LifeV::Operators::FSIApplyOperator> M_applyOperatorResidual;

	// Operator to apply Jacobian of the FSI in the nonconforming case
	boost::shared_ptr<LifeV::Operators::FSIApplyOperatorNonConforming> M_applyOperatorJacobianNonConforming;

	// Preconditioner operator
	boost::shared_ptr<LifeV::Operators::DirichletNeumannPreconditioner> M_prec;

	// Epetra Operator needed to solve the linear system
	boost::shared_ptr<Operators::InvertibleOperator> M_invOper;

	// Parameter list solver
	parameterListPtr_Type M_pListLinSolver;

	// output txt files
	UInt M_NewtonIter;
	bool M_useShapeDerivatives;
	bool M_printResiduals;
	bool M_printSteps;
	std::ofstream M_outputTimeStep;
	std::ofstream M_outputResiduals;
	std::ofstream M_outputSteps;
    std::ofstream M_outputLinearIterations;
    std::ofstream M_outputPreconditionerComputation;
    std::ofstream M_outputTimeLinearSolver;
    
	// Extrapolation of the initial guess for Newton
	bool M_extrapolateInitialGuess;
	UInt M_orderExtrapolationInitialGuess;
	boost::shared_ptr<TimeAndExtrapolationHandler> M_extrapolationSolution;

	// paritioned meshes
	bool M_usePartitionedMeshes;
	boost::shared_ptr<std::map<UInt, UInt> > M_localDofMap;

	bool M_subiterateFluidDirichlet;

	bool M_considerGravity;
	Real M_gravity;
	Real M_gravityDirection;

	bool M_moveMesh;

	bool M_nonconforming;
	interpolationPtr_Type M_FluidToStructureInterpolant;
	interpolationPtr_Type M_StructureToFluidInterpolant;

	matrixPtr_Type M_interface_mass_structure;
	matrixPtr_Type M_interface_mass_structure_robin;
	matrixPtr_Type M_interface_mass_fluid;

	bool M_lambda_num_structure;

	mapPtr_Type M_lagrangeMapScalar;
	vectorPtr_Type M_numerationInterfaceFluid;
	vectorPtr_Type M_numerationInterfaceStructure;

	bool M_useMasses;
    
    // Members for the restarter
    
    bool M_restart;
    
    boost::shared_ptr<ExporterHDF5<mesh_Type > > M_importerFluid;     // Import just from solution in hdf5 format! I do not know
    boost::shared_ptr<ExporterHDF5<mesh_Type > > M_importerStructure; // if the import with other filters works.

	// To handle the post-processing
	int M_saveEvery;
	int M_counterSaveEvery;

	vectorPtr_Type M_rhsStructureVelocity;

	basePrecPtr_Type M_precPtr;
	bool 			 M_precPtrBuilt;
	bool 			 M_linearElasticity;

	bool             M_disregardRestart; // if true disregards correct exporter as function of BDF
	bool             M_prescribeInflowFlowrate;
	vectorPtr_Type   M_dsk;
	bool             M_useBDF;
	boost::shared_ptr<BDFSecondOrderDerivative> M_structureTimeAdvanceBDF;
	UInt             M_orderBDFSolid;
};

} // end namespace LifeV

#endif // end SETTINGS_H
