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
#include <lifev/navier_stokes_blocks/solver/NavierStokesSolver.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialLinear.hpp>

#include <lifev/fsi_blocks/solver/ALESolver.hpp>

// Expression template FE space
#include <lifev/eta/fem/ETFESpace.hpp>

// time advance for the structure
#include <lifev/core/fem/TimeAdvance.hpp>
#include <lifev/core/fem/TimeAdvanceNewmark.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include <lifev/core/fem/TimeAndExtrapolationHandler.hpp>

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

// Includes needed by interpolation when using nonconforming meshes
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <lifev/core/interpolation/RBFInterpolation.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>

#include <lifev/fsi_blocks/solver/FSIApplyOperatorNonConforming.hpp>

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

	typedef MatrixEpetra<Real> matrix_Type;
	typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

	typedef Teuchos::ParameterList parameterList_Type;
	typedef boost::shared_ptr<parameterList_Type> parameterListPtr_Type;

	typedef RBFInterpolation<mesh_Type>           interpolation_Type;
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

    void setBoundaryConditions ( const bcPtr_Type& fluidBC, const bcPtr_Type& fluidBC_residual, const bcPtr_Type& structureBC, const bcPtr_Type& aleBC);

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

    void setBoundaryConditionsPCD ( const bcPtr_Type& pcdBC);

    void setGravity ( const Real& gravity, const Real& gravity_direction);

//@}

private:

    void moveMesh ( const VectorEpetra& displacement );

    void buildMonolithicMap ( );

    void createStructureFESpaces ( );

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

    void getRhsStructure ( );

    void applyBCstructure ( );

    void applyBCsolution(vectorPtr_Type& M_solution);

    void applyBCresidual(VectorEpetra& residual);

    void setSolversOptions(const Teuchos::ParameterList& solversOptions);

    void updateRhsCouplingVelocities ( );

    void updateRhsCouplingVelocities_nonconforming ( );

    void structureToInterface (vector_Type& VectorOnGamma, const vector_Type& VectorOnStructure);

    void initializeExtrapolation( );

    void assembleStructureInterfaceMass ( );

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
    boost::shared_ptr<NavierStokesSolver> M_fluid;
    boost::shared_ptr<StructuralOperator<mesh_Type> > M_structure;
    boost::shared_ptr<StructuralConstitutiveLawData> M_dataStructure;
    boost::shared_ptr<ALESolver> M_ale;

    // time advance for the structure
    timeAdvancePtr_Type M_structureTimeAdvance;
    timeAdvancePtr_Type M_aleTimeAdvance;
    boost::shared_ptr<TimeAndExtrapolationHandler> M_fluidTimeAdvance;

    // boundary conditions
    bcPtr_Type M_fluidBC;
    bcPtr_Type M_fluidBC_residual;
    bcPtr_Type M_structureBC;
    bcPtr_Type M_aleBC;
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
	vectorPtr_Type M_structureDisplacement;

	// Monolithic map
	boost::shared_ptr<map_Type> M_monolithicMap;

	vectorPtr_Type M_rhsFluid;
	vectorPtr_Type M_rhsStructure;
	vectorPtr_Type M_rhsCouplingVelocities;
	vectorPtr_Type M_solution;

	vectorPtr_Type M_u_star;
	vectorPtr_Type M_w_star;
	vectorPtr_Type M_beta_star;
	vectorPtr_Type M_rhs_velocity;

	matrixPtr_Type M_matrixStructure;

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

	// Extrapolation of the initial guess for Newton
	bool M_extrapolateInitialGuess;
	UInt M_orderExtrapolationInitialGuess;
	boost::shared_ptr<TimeAndExtrapolationHandler> M_extrapolationSolution;

	// paritioned meshes
	bool M_usePartitionedMeshes;
	boost::shared_ptr<std::map<UInt, UInt> > M_localDofMap;

	//! BCs for the PCD block Fp
	bcPtr_Type M_pcdBC;
	bool M_subiterateFluidDirichlet;

	bool M_considerGravity;
	Real M_gravity;
	Real M_gravityDirection;

	bool M_moveMesh;

	bool M_nonconforming;
	interpolationPtr_Type M_FluidToStructureInterpolant;
	interpolationPtr_Type M_StructureToFluidInterpolant;

	matrixPtr_Type M_interface_mass_structure;
	matrixPtr_Type M_interface_mass_fluid;
};

} // end namespace LifeV

#endif // end SETTINGS_H
