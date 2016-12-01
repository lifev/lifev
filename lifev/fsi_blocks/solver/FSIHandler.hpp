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

    void setDatafile(const GetPot& dataFile);

    void readMeshes( );

    void partitionMeshes ( );

    void readPartitionedMeshes( );

    void setup ( );

    void setBoundaryConditions ( const bcPtr_Type& fluidBC, const bcPtr_Type& structureBC, const bcPtr_Type& aleBC );

    void setFluidInterfaceBoundaryConditions ( const bcPtr_Type& interfaceFluidBC ) { M_interfaceFluidBC = interfaceFluidBC; };

    // update all the bc handlers
    void updateBoundaryConditions( );

    void initializeTimeAdvance ( );

    void buildInterfaceMaps ( );

    void assembleCoupling ( );

    void solveFSIproblem();

    void evalResidual(vector_Type& residual, const vector_Type& solution, const UInt iter_newton);

    void solveJac( vector_Type& increment, const vector_Type& residual, const Real linearRelTol );

    void setParameterLists( );

    void setGravity ( const Real& gravity, const Real& gravity_direction);

    Real getTime () { return M_time; };

    Real getStartTime () { return M_t_zero; };

    void setTime (const Real& time ) { M_time = time; };

    Real getTimeStep () { return M_dt; };

    Real getEndTime () { return M_t_end; };

    boost::shared_ptr<NavierStokesSolverBlocks> getFluid() { return M_fluid; };

    void intializeTimeLoop ( );

    void solveTimeStep( );

    void postprocessResults( const int& time_step_count );

    void finalizeExporters ( );

    vectorPtr_Type getFluidVelocity() { return M_fluidVelocity; };

    vectorPtr_Type getFluidPressure() { return M_fluidPressure; };

//@}

private:

    void moveMesh ( const VectorEpetra& displacement );

    void buildMonolithicMap ( );

    void setupStructure ( );

    void createAleFESpace();

    // update the bc handler
    void updateBCHandler( bcPtr_Type & bc );

    void createInterfaceMaps ( std::map<ID, ID> const& locDofMap );

    void constructInterfaceMap ( const std::map<ID, ID>& locDofMap, const UInt subdomainMaxId);

    void setupExporters( );

    void instantiateExporter( boost::shared_ptr< Exporter<mesh_Type > > & exporter,
			  	  	  	  	  const meshPtr_Type& localMesh,
			  	  	  	  	  const std::string& outputFileName );

    void updateSystem ( );

    void initializeApplyOperatorResidual ( );

    void initializeApplyOperatorJacobian ( );

    void getMatrixStructure ( );

    void get_structure_coupling_velocities ( );

    void applyBCstructure ( );

    void applyBCsolution(vectorPtr_Type& M_solution);

    void applyBCresidual(VectorEpetra& r_u, VectorEpetra& r_ds, VectorEpetra& r_df);

    void setSolversOptions(const Teuchos::ParameterList& solversOptions);

    void updateRhsCouplingVelocities ( );

    void updateRhsCouplingVelocities_nonconforming ( );

    void structureToInterface (vector_Type& VectorOnGamma, const vector_Type& VectorOnStructure);

    void initializeExtrapolation( );

    void assembleStructureInterfaceMass ( );

    void applyInverseFluidMassOnGamma ( const vectorPtr_Type& lambda, vectorPtr_Type& strongLambda );

    void solveFSIproblemAorta();

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
