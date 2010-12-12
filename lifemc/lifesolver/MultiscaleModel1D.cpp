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
 *  @file
 *  @brief File containing the MultiScale Model 1D
 *
 *  @version 1.1
 *  @date 26-02-2010
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *
 *  @version 1.2 and subsequents
 *  @date 23-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include "MultiscaleModel1D.hpp"

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Model_1D::MS_Model_1D() :
        MS_PhysicalModel                          (),
#ifdef HAVE_HDF5
        M_exporter                     ( new IOFile_Type() ),
        M_importer                     ( new IOFile_Type() ),
        M_exporterMesh                 ( new Mesh_Type() ),
#endif
#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
        M_linearBC                     ( new BC_Type() ),
        M_linearSolution               ( new Solution_Type() ),
        M_bcPreviousTimeSteps          (),
        M_bcBaseDelta                  (),
        M_bcDelta                      (),
        M_bcDeltaType                  (),
        M_bcDeltaSide                  (),
#endif
        M_data                         ( new Data_Type() ),
        M_bc                           ( new BCInterface_Type() ),
        M_physics                      (),
        M_flux                         (),
        M_source                       (),
        M_solver                       ( new Solver_Type() ),
        M_linearSolver                 (),
        M_FESpace                      (),
        M_solution_tn                  ( new Solution_Type() ),
        M_solution                     ( new Solution_Type() )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::MS_Model_1D() \n";
#endif

    M_type = OneDimensional;

    //Define the maps of the OneDimensionalModel objects
    OneDimensionalModel_MapsDefinition();

    //Register the objects
    factoryOneDimensionalPhysics_Type::instance().registerProduct( OneD_LinearPhysics,    &createOneDimensionalPhysicsLinear );
    factoryOneDimensionalPhysics_Type::instance().registerProduct( OneD_NonLinearPhysics, &createOneDimensionalPhysicsNonLinear );

    factoryOneDimensionalFlux_Type::instance().registerProduct(    OneD_LinearFlux,       &createOneDimensionalFluxLinear );
    factoryOneDimensionalFlux_Type::instance().registerProduct(    OneD_NonLinearFlux,    &createOneDimensionalFluxNonLinear );

    factoryOneDimensionalSource_Type::instance().registerProduct(  OneD_LinearSource,     &createOneDimensionalSourceLinear );
    factoryOneDimensionalSource_Type::instance().registerProduct(  OneD_NonLinearSource,  &createOneDimensionalSourceNonLinear );
}

// ===================================================
// MultiScale PhysicalModel Virtual Methods
// ===================================================
void
MS_Model_1D::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupData( ) \n";
#endif

    // Preliminary setup of the communicator
#ifdef EPETRA_MPI
    MPI_Comm LocalComm;
    MPI_Comm_split( ( dynamic_cast<Epetra_MpiComm*> ( &(*M_comm) ) )->Comm(), M_comm->MyPID(), M_comm->MyPID(), &LocalComm );
    M_comm.reset( new Epetra_MpiComm( LocalComm ) );
#else
    M_comm.reset( new Epetra_SerialComm() );
#endif

    MS_PhysicalModel::setupData( fileName );

    GetPot dataFile( fileName );

    M_data->setup( dataFile );
    if ( M_globalData.get() )
        setupGlobalData( fileName );

    //1D Model Physics
    M_physics = Physics_PtrType( factoryOneDimensionalPhysics_Type::instance().createObject( M_data->physicsType() ) );
    M_physics->setData( M_data );

    //1D Model Flux
    M_flux = Flux_PtrType( factoryOneDimensionalFlux_Type::instance().createObject( M_data->fluxType() ) );
    M_flux->setPhysics( M_physics );

    //1D Model Source
    M_source = Source_PtrType( factoryOneDimensionalSource_Type::instance().createObject( M_data->sourceType() ) );
    M_source->setPhysics( M_physics );

    //Linear Solver
    M_linearSolver.reset( new LinearSolver_Type( M_comm ) );
    M_linearSolver->setUpPrec        ( dataFile, "1D_Model/prec" );
    M_linearSolver->setDataFromGetPot( dataFile, "1D_Model/solver" );
    M_linearSolver->setParameters();

    //1D Model Solver
    M_solver->setCommunicator( M_comm );
    M_solver->setProblem( M_physics, M_flux, M_source );
    M_solver->setLinearSolver( M_linearSolver );

    //BC - We need to create the BCHandler before using it
    M_bc->createHandler();
    //M_bc->fillHandler( fileName, "1D_Model" );

    //Exporters
    M_data->setPostprocessingDirectory( MS_ProblemFolder );
    M_data->setPostprocessingFile( "Step_" + number2string( MS_ProblemStep ) + "_Model_" + number2string( M_ID ) );

#ifdef HAVE_HDF5
    M_exporter->setDataFromGetPot( dataFile );
    M_exporter->setPrefix( "Step_" + number2string( MS_ProblemStep ) + "_Model_" + number2string( M_ID ) );
    M_exporter->setDirectory( MS_ProblemFolder );

    M_exporterMesh->setup( M_data->length(), M_data->numberOfElements() );

    M_importer->setDataFromGetPot( dataFile );
    M_importer->setPrefix( "Step_" + number2string( MS_ProblemStep - 1 ) + "_Model_" + number2string( M_ID ) );
    M_importer->setDirectory( MS_ProblemFolder );
#endif

}

void
MS_Model_1D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupProblem() \n";
#endif

    //FEspace
    setupFESpace();

    //Setup solution
    M_solver->setupSolution( *M_solution );
    M_solver->setupSolution( *M_solution_tn );

    //Set default BC (has to be called after setting other BC)
    M_bc->handler()->setDefaultBC();
    M_bc->setPhysicalSolver( M_solver );
    M_bc->setSolution( M_solution );
    M_bc->setFluxSource( M_flux, M_source );

#ifdef HAVE_HDF5
    //Post-processing
    M_exporter->setMeshProcId( M_exporterMesh, M_comm->MyPID() );

    //M_exporter->addVariable( ExporterData::Scalar, "Solid Area",      (*M_Solution)["A"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_exporter->addVariable( ExporterData::Scalar, "Area ratio",      (*M_solution)["A/A0-1"], static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_exporter->addVariable( ExporterData::Scalar, "Fluid Flow Rate", (*M_solution)["Q"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    //M_exporter->addVariable( ExporterData::Scalar, "W1",              (*M_Solution)["W1"],   static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    //M_exporter->addVariable( ExporterData::Scalar, "W2",              (*M_Solution)["W2"],   static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_exporter->addVariable( ExporterData::Scalar, "Fluid Pressure",  (*M_solution)["P"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
#endif

#ifdef HAVE_MATLAB_POSTPROCESSING
    //Matlab post-processing
    M_solver->resetOutput( *M_solution );
#endif

    //Setup solution
    initializeSolution();

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
    {
        createLinearBC();
        updateLinearBC( *M_solution );
        setupLinearModel();
    }
#endif

}

void
MS_Model_1D::buildSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::BuildSystem() \n";
#endif

    //M_Data->showMe();
    M_solver->buildConstantMatrices();

    // Update previous solution
    updateSolution( *M_solution, *M_solution_tn );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
        updateLinearModel();
#endif

}

void
MS_Model_1D::updateSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::UpdateSystem() \n";
#endif

    // Update previous solution
    updateSolution( *M_solution, *M_solution_tn );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
        updateLinearModel();
#endif

}

void
MS_Model_1D::solveSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SolveSystem() \n";
#endif

    solve( *M_bc->handler(), *M_solution );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
        updateLinearBC( *M_solution );
#endif

}

void
MS_Model_1D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SaveSolution() \n";
#endif

#ifdef HAVE_HDF5
    M_exporter->postProcess( M_data->dataTime()->getTime() );

    if ( M_data->dataTime()->isLastTimeStep() )
        M_exporter->CloseFile();
#endif

#ifdef HAVE_MATLAB_POSTPROCESSING
    //Matlab post-processing
    M_solver->postProcess( *M_solution );
#endif

}

void
MS_Model_1D::showMe()
{
    if ( M_displayer->isLeader() )
    {
        MS_PhysicalModel::showMe();

        std::cout << "FE order            = " << "P1" << std::endl
                  << "DOF                 = " << M_data->mesh()->numPoints() << std::endl << std::endl;

        std::cout << "maxH                = " << M_data->mesh()->maxH() << std::endl
                  << "meanH               = " << M_data->mesh()->meanH() << std::endl << std::endl;
    }
}

// ===================================================
// Methods
// ===================================================
#ifdef JACOBIAN_WITH_FINITEDIFFERENCE

void
MS_Model_1D::setupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupLinearModel( ) \n";
#endif

    // Define BCFunction for linear problem
    M_bcBaseDelta.setFunction( boost::bind( &MS_Model_1D::bcFunctionDelta, this, _1 ) );

    // The linear BCHandler is a copy of the original BCHandler with the LinearSolution instead of the true solution
    //M_LinearBC.reset( new bc_Type( *M_bc->handler() ) ); // COPY CONSTRUCTOR NOT WORKING

    //Set left and right BC + default BC
    M_linearBC->setBC( OneD_left, OneD_first, M_bc->handler()->BC( OneD_left )->type( OneD_first ),
                       M_bc->handler()->BC( OneD_left )->BCFunction( OneD_first ) );

    M_linearBC->setBC( OneD_right, OneD_first, M_bc->handler()->BC( OneD_right )->type( OneD_first ),
                       M_bc->handler()->BC( OneD_right )->BCFunction( OneD_first ) );

    M_linearBC->setDefaultBC();

    // Solution for the linear problem (this does not change anything in the solver)
    M_solver->setupSolution( *M_linearSolution );
    M_linearBC->setSolution( M_linearSolution );
    M_linearBC->setFluxSource( M_flux, M_source );
}

void
MS_Model_1D::updateLinearModel()
{
    // The couplings should use the same value for the time interpolation order
    UInt timeInterpolationOrder( std::max( M_couplings[0]->timeInterpolationOrder(), M_couplings[1]->timeInterpolationOrder() ) );

    UInt containerSize( M_bcPreviousTimeSteps.size() );

    // If we have not yet enough samples for interpolation, we add a new one
    if ( containerSize <= timeInterpolationOrder )
    {
        ++containerSize;
        M_bcPreviousTimeSteps.push_back( M_bcPreviousTimeSteps[0] );
    }

    // Updating database
    for ( UInt i(1) ; i < containerSize ; ++i )
        M_bcPreviousTimeSteps[containerSize-i] = M_bcPreviousTimeSteps[containerSize-i-1];
}

void
MS_Model_1D::solveLinearModel( bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SolveLinearModel() \n";
#endif

    if ( !solveLinearSystem )
        return;

    imposePerturbation();

    solve( *M_linearBC, *M_linearSolution, "L1D-" );

    resetPerturbation();

    //This flag avoid recomputation of the same system
    solveLinearSystem = false;
}

#endif

// ===================================================
// Get Methods (couplings)
// ===================================================
MS_Model_1D::BCInterface_Type&
MS_Model_1D::bcInterface() const
{
    return *M_bc;
}

Real
MS_Model_1D::boundaryDensity( const BCFlag& /*flag*/ ) const
{
    return M_data->densityRho();
}

Real
MS_Model_1D::boundaryViscosity( const BCFlag& /*flag*/ ) const
{
    return M_data->viscosity();
}

Real
MS_Model_1D::boundaryArea( const BCFlag& flag ) const
{
    return M_solver->boundaryValue( *M_solution, OneD_A, flagConverter( flag ) );
}

Real
MS_Model_1D::boundaryFlowRate( const BCFlag& flag ) const
{
    return M_solver->boundaryValue( *M_solution, OneD_Q, flagConverter( flag ) );
}

Real
MS_Model_1D::boundaryPressure( const BCFlag& flag ) const
{
    return M_solver->boundaryValue( *M_solution, OneD_P, flagConverter( flag ) );
}

Real
MS_Model_1D::boundaryDynamicPressure( const BCFlag& flag ) const
{
    return 0.5 * boundaryDensity( flag ) * std::pow( boundaryFlowRate( flag ) / boundaryArea( flag ), 2 );
}

Real
MS_Model_1D::boundaryStress( const BCFlag& flag, const stress_Type& stressType ) const
{
    switch ( stressType )
    {
    case StaticPressure:
    {
        return boundaryPressure( flag );
    }

    case TotalPressure:
    {
        return boundaryPressure( flag ) + boundaryDynamicPressure( flag ) * ( ( boundaryFlowRate( flag ) > 0.0 ) ? 1 : -1 );
    }

    default:
    {
        std::cout << "ERROR: Invalid stress type [" << Enum2String( stressType, MS_stressesMap ) << "]" << std::endl;

        return 0.0;
    }
    }
}

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE

Real
MS_Model_1D::boundaryDeltaFlowRate( const BCFlag& flag, bool& solveLinearSystem )
{
    OneD_BCSide bcSide = flagConverter( flag );

    solveLinearModel( solveLinearSystem );

    Real Q      = M_solver->boundaryValue( *M_solution, OneD_Q, bcSide );
    Real Qdelta = M_solver->boundaryValue( *M_linearSolution, OneD_Q, bcSide );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::GetBoundaryDeltaFlowRate( flag, solveLinearSystem ) \n";
    Debug( 8130 ) << "Q:          " << Q << "\n";
    Debug( 8130 ) << "Qdelta:     " << Qdelta << "\n";
#endif

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA

    if ( M_bcDeltaType == OneD_A )
    {
        // dQ/dP
        return ( (Qdelta - Q) / M_bcDelta ) * M_physics->dAdP( M_solver->boundaryValue( *M_solution, OneD_P, bcSide ), 0 );
    }
    else
    {
        // dQ/dQ
        return (Qdelta - Q) / M_bcDelta;
    }

#else

    return (Qdelta - Q) / M_bcDelta;

#endif

}

Real
MS_Model_1D::boundaryDeltaPressure( const BCFlag& flag, bool& solveLinearSystem )
{
    OneD_BCSide bcSide = flagConverter( flag );

    solveLinearModel( solveLinearSystem );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA

    Real A      = M_solver->boundaryValue( *M_solution, OneD_A, bcSide );
    Real Adelta = M_solver->boundaryValue( *M_linearSolution, OneD_A, bcSide );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::GetBoundaryDeltaPressure( flag, solveLinearSystem ) \n";
    Debug( 8130 ) << "A:          " << A <<  "\n";
    Debug( 8130 ) << "Adelta:     " << Adelta <<  "\n";
#endif

    if ( M_bcDeltaType == OneD_A )
    {
        // dP/dP
        return ( (Adelta - A) / M_bcDelta ) * M_physics->dPdA( M_solver->boundaryValue( *M_solution, OneD_A, bcSide ), 0 )
               * M_physics->dAdP( M_solver->boundaryValue( *M_solution, OneD_P, bcSide ), 0 );
    }
    else
    {
        // dP/dQ
        return ( (Adelta - A) / M_bcDelta ) * M_physics->dPdA( M_solver->boundaryValue( *M_solution, OneD_A, bcSide ), 0 );
    }

#else

    Real P      = M_solver->boundaryValue( *M_solution, OneD_P, bcSide );
    Real Pdelta = M_solver->boundaryValue( *M_linearSolution, OneD_P, bcSide );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::GetBoundaryDeltaPressure( flag, solveLinearSystem ) \n";
    Debug( 8130 ) << "P:          " << P <<  "\n";
    Debug( 8130 ) << "Pdelta:     " << Pdelta <<  "\n";
#endif

    return (Pdelta - P) / M_bcDelta;

#endif

}

#else

Real
MS_Model_1D::boundaryDeltaFlowRate( const BCFlag& flag, bool& /*solveLinearSystem*/ )
{
    return tangentProblem( flagConverter( flag ), OneD_Q );
}

Real
MS_Model_1D::boundaryDeltaPressure( const BCFlag& flag, bool& /*solveLinearSystem*/ )
{
    return tangentProblem( flagConverter( flag ), OneD_P );
}

#endif

Real
MS_Model_1D::boundaryDeltaDynamicPressure( const BCFlag& flag, bool& solveLinearSystem )
{
    // Untested
    return boundaryDensity( flag ) * boundaryDeltaFlowRate( flag, solveLinearSystem ) * boundaryFlowRate( flag ) / ( boundaryArea( flag ) * boundaryArea( flag ) );
}

Real
MS_Model_1D::boundaryDeltaStress( const BCFlag& flag, bool& solveLinearSystem, const stress_Type& stressType )
{
    switch ( stressType )
    {
    case StaticPressure:
    {
        return boundaryDeltaPressure( flag, solveLinearSystem );
    }

    case TotalPressure:
    {
        return boundaryDeltaPressure( flag, solveLinearSystem ) + boundaryDeltaDynamicPressure( flag, solveLinearSystem ); //Verify the sign of DynamicPressure contribute!
    }

    default:
    {
        std::cout << "ERROR: Invalid stress type [" << Enum2String( stressType, MS_stressesMap ) << "]" << std::endl;

        return 0.0;
    }
    }
}

// ===================================================
// Get Methods
// ===================================================
MS_Model_1D::BC_Type&
MS_Model_1D::bc() const
{
    return *(M_bc->handler());
}


MS_Model_1D::Data_Type&
MS_Model_1D::data() const
{
    return *M_data;
}

MS_Model_1D::Physics_PtrType
MS_Model_1D::physics() const
{
    return M_physics;
}

MS_Model_1D::Flux_PtrType
MS_Model_1D::flux() const
{
    return M_flux;
}

MS_Model_1D::Source_PtrType
MS_Model_1D::source() const
{
    return M_source;
}

MS_Model_1D::FESpace_PtrType
MS_Model_1D::FESpace() const
{
    return M_FESpace;
}

MS_Model_1D::Solver_PtrType
MS_Model_1D::solver() const
{
    return M_solver;
}

const MS_Model_1D::Solution_PtrType&
MS_Model_1D::solution() const
{
    return M_solution;
}

const MS_Model_1D::Vector_PtrType&
MS_Model_1D::solution( const std::string& quantity) const
{
    return (*M_solution)[quantity];
}

// ===================================================
// Private Methods
// ===================================================
void
MS_Model_1D::setupGlobalData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupGlobalData( fileName ) \n";
#endif

    GetPot dataFile( fileName );

    //Global data time
    M_data->setDataTime( M_globalData->dataTime() );

    //Global physical quantities
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/density" ) )
        M_data->setDensity( M_globalData->fluidDensity() );
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/viscosity" ) )
        M_data->setViscosity( M_globalData->fluidViscosity() );
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/externalPressure" ) )
        M_data->setExternalPressure( M_globalData->fluidReferencePressure() );
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/densityWall" ) )
        M_data->setDensityWall( M_globalData->structureDensity() );
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/poisson" ) )
        M_data->setPoisson( M_globalData->structurePoissonCoefficient() );
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/young" ) )
        M_data->setYoung( M_globalData->structureYoungModulus() );

    //After changing some parameters we need to update the coefficients
    M_data->updateCoefficients();
}

void
MS_Model_1D::setupFESpace()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupFEspace() \n";
#endif

    //Transform mesh
    boost::array< Real, NDIM > NullTransformation;
    NullTransformation[0] = 0.;
    NullTransformation[1] = 0.;
    NullTransformation[2] = 0.;

    //The real mesh can be only scaled due to OneDimensionalModel_Solver conventions
    M_data->mesh()->transformMesh( M_geometryScale, NullTransformation, NullTransformation ); // Scale the x dimension

    for ( UInt i(0); i < M_data->numberOfNodes() ; ++i )
        M_data->setArea0( M_data->area0( i ) * M_geometryScale[1] * M_geometryScale[2], i );  // Scale the area (y-z dimensions)

    //After changing some parameters we need to update the coefficients
    M_data->updateCoefficients();

#ifdef HAVE_HDF5
    //The mesh for the post-processing can be rotated
    M_exporterMesh->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );
#endif

    //Setup FESpace
    const RefFE*    refFE = &feSegP1;
    const QuadRule* qR    = &quadRuleSeg3pt;
    const QuadRule* bdQr  = &quadRuleSeg1pt;

//    const RefFE*    refFE = &feSegP2;
//    const QuadRule* qR    = &quadRuleSeg3pt;
//    const QuadRule* bdQr  = &quadRuleSeg1pt;

    M_FESpace.reset( new FESpace_Type( M_data->mesh(), *refFE, *qR, *bdQr, 1, M_comm ) );
    M_solver->setFESpace( M_FESpace );
}

void
MS_Model_1D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::InitializeSolution() \n";
#endif

    if ( MS_ProblemStep > 0 )
    {
        M_importer->setMeshProcId( M_exporterMesh, M_comm->MyPID() );

        //M_exporter->addVariable( ExporterData::Scalar, "Solid Area",      (*M_Solution)["A"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        M_importer->addVariable( ExporterData::Scalar, "Area ratio",      (*M_solution)["A/A0-1"], static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        M_importer->addVariable( ExporterData::Scalar, "Fluid Flow Rate", (*M_solution)["Q"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        //M_importer->addVariable( ExporterData::Scalar, "W1",              (*M_Solution)["W1"],   static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        //M_importer->addVariable( ExporterData::Scalar, "W2",              (*M_Solution)["W2"],   static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        M_importer->addVariable( ExporterData::Scalar, "Fluid Pressure",  (*M_solution)["P"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );

        // Import
        M_exporter->setStartIndex( M_importer->importFromTime( M_data->dataTime()->getInitialTime() ) + 1 );

        // Compute A from AreaRatio
        M_solver->computeArea( *M_solution );

        // Compute W1 and W2 from A and Q
        M_solver->computeW1W2( *M_solution );
    }
    else
        M_solver->initialize( *M_solution );
}

void
MS_Model_1D::updateSolution( const Solution_Type& solution1, Solution_Type& solution2 )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::UpdateSolution( solution1, solution2 ) \n";
#endif

    // Here we make a true copy (not a copy of the shared_ptr, but a copy of its content)
    for ( Solution_ConstIterator i = solution1.begin() ; i != solution1.end() ; ++i )
        *solution2[i->first] = *i->second;
}

void
MS_Model_1D::solve( BC_Type& bc, Solution_Type& solution, const std::string& solverType )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::Solve() \n";
#endif

    // Re-initialize solution
    updateSolution( *M_solution_tn, solution );

    // Subiterate to respect CFL
    UInt SubiterationNumber(1);
    Real timeStep = M_data->dataTime()->getTimeStep();

    Real CFL = M_solver->computeCFL( solution, M_data->dataTime()->getTimeStep() );
    if ( CFL > M_data->CFLmax() )
    {
        SubiterationNumber = std::ceil( CFL / M_data->CFLmax() );
        timeStep /= SubiterationNumber;
    }

    if ( M_displayer->isLeader() )
        std::cout << solverType << "  CFL                                      " << CFL*timeStep/M_data->dataTime()->getTimeStep() << std::endl;

    for ( UInt i(1) ; i <= SubiterationNumber ; ++i )
    {
        if ( M_displayer->isLeader() )
        {
            std::cout << solverType << "  Subiteration                             " << i << "/" << SubiterationNumber << std::endl;
            std::cout << solverType << "  Time                                     " <<  M_data->dataTime()->getPreviousTime() + i*timeStep << std::endl;
        }
        //bc.updateOperatorVariables();

        M_solver->updateRHS( solution, timeStep );
        M_solver->iterate( bc, solution, M_data->dataTime()->getPreviousTime() + i*timeStep, timeStep );
    }
}

OneD_BCSide
MS_Model_1D::flagConverter( const BCFlag& flag ) const
{
    return (flag == 0) ? OneD_left : OneD_right;
}

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE

void
MS_Model_1D::createLinearBC()
{
    // Allocating the correct space
    M_bcPreviousTimeSteps.reserve( std::max( M_couplings[0]->timeInterpolationOrder(), M_couplings[1]->timeInterpolationOrder() ) );

    // Create bcSide map
    std::map< OneD_BCSide, std::map< OneD_BC, Real > > bcSideMap;
    M_bcPreviousTimeSteps.push_back( bcSideMap );

    // Create bcType map
    std::map< OneD_BC, Real > bcTypeMap;
    M_bcPreviousTimeSteps[0][OneD_left]  = bcTypeMap;
    M_bcPreviousTimeSteps[0][OneD_right] = bcTypeMap;
}

void
MS_Model_1D::updateLinearBC( const Solution_Type& solution )
{
    M_bcPreviousTimeSteps[0][OneD_left][OneD_A]  = M_solver->boundaryValue( solution, OneD_A, OneD_left );
    M_bcPreviousTimeSteps[0][OneD_left][OneD_P]  = M_solver->boundaryValue( solution, OneD_P, OneD_left );
    M_bcPreviousTimeSteps[0][OneD_left][OneD_Q]  = M_solver->boundaryValue( solution, OneD_Q, OneD_left );
    M_bcPreviousTimeSteps[0][OneD_right][OneD_A] = M_solver->boundaryValue( solution, OneD_A, OneD_right );
    M_bcPreviousTimeSteps[0][OneD_right][OneD_P] = M_solver->boundaryValue( solution, OneD_P, OneD_right );
    M_bcPreviousTimeSteps[0][OneD_right][OneD_Q] = M_solver->boundaryValue( solution, OneD_Q, OneD_right );
}

void
MS_Model_1D::imposePerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::Perturbation() \n";
#endif

    for ( MS_CouplingsVector_ConstIterator i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            // Find the side to perturb and apply the perturbation
            M_bcDeltaSide = flagConverter( ( *i )->flag( ( *i )->modelGlobalToLocalID( M_ID ) ) );
            M_linearBC->BC( M_bcDeltaSide )->setBCFunction( OneD_first, M_bcBaseDelta );

            // Compute the range
            M_bcDeltaType = M_linearBC->BC( M_bcDeltaSide )->type( OneD_first );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA
            // We replace pressure BC with area BC for the perturbed problem
            if ( M_bcDeltaType == OneD_P )
            {
                M_linearBC->BC( M_bcDeltaSide )->setType( OneD_first, OneD_A );
                M_bcDeltaType = OneD_A;
            }
#endif

            //M_BCDelta = ( *i )->residual()[( *i )->perturbedCoupling()] * 10000;
            //M_BCDelta = ( M_BCDelta[1] - M_BCDelta[0] ) / 100;

            //if ( std::abs( M_BCDelta ) < 1e-6 || std::abs( M_BCDelta ) > 1e6 )
            switch ( M_bcDeltaType )
            {
            case OneD_A:

                M_bcDelta = M_data->jacobianPerturbationArea();
                break;

            case OneD_Q:

                M_bcDelta = M_data->jacobianPerturbationFlowRate();

                break;

            case OneD_P:

                M_bcDelta = 5; //M_Data->jacobianPerturbationPressure();

                break;

            default:
                std::cout << "Warning: bcType \"" << M_bcDeltaType << "\"not available!" << std::endl;
            }

            break;
        }

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "BCDelta:    " << M_bcDelta << "\n";
#endif

}

void
MS_Model_1D::resetPerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::ResetPerturbation() \n";
#endif

    M_linearBC->BC( M_bcDeltaSide )->setBCFunction( OneD_first, M_bc->handler()->BC( M_bcDeltaSide )->BCFunction( OneD_first ) );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA
    // Restoring the original BC
    if ( M_bcDeltaType == OneD_A )
        M_linearBC->BC( M_bcDeltaSide )->setType( OneD_first, OneD_P );
#endif

}

Real
MS_Model_1D::bcFunctionDelta( const Real& t )
{
    // Lagrange interpolation
    Real bcValue(0);
    Real base(1);

    // Time container for interpolation
    std::vector< Real > timeContainer( M_bcPreviousTimeSteps.size(), 0 );
    for ( UInt i(0) ; i < M_bcPreviousTimeSteps.size() ; ++i )
        timeContainer[i] = M_globalData->dataTime()->getTime() - i * M_globalData->dataTime()->getTimeStep();

    for ( UInt i(0) ; i < M_bcPreviousTimeSteps.size() ; ++i )
    {
        base = 1;
        for ( UInt j(0) ; j < M_bcPreviousTimeSteps.size() ; ++j )
            if ( j != i )
                base *= (t - timeContainer[j]) / (timeContainer[i] - timeContainer[j]);

        if ( i == 0 )
            bcValue += ( M_bcDelta + M_bcPreviousTimeSteps[i][M_bcDeltaSide][M_bcDeltaType] ) * base;
        else
            bcValue += M_bcPreviousTimeSteps[i][M_bcDeltaSide][M_bcDeltaType] * base;
    }

    return bcValue;
}

#else

Real
MS_Model_1D::tangentProblem( const OneD_BCSide& bcOutputSide, const OneD_BC& bcOutputType )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::TangentProblem( bcOutputSide, bcOutputType ) \n";
#endif

    Real jacobianCoefficient(0);

    for ( MS_CouplingsVector_ConstIterator i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            // Find the perturbed side
            OneD_BCSide bcSide = flagConverter( ( *i )->flag( ( *i )->modelGlobalToLocalID( M_ID ) ) );

            // Perturbation has no effect on the other sides (which also means that dQ/dQ and dP/dP are always zero)
            if ( bcSide != bcOutputSide )
                break;

            // Compute the eigenvectors
            Container2D_Type eigenvalues, leftEigenvector1, leftEigenvector2;
            M_solver->BoundaryEigenValuesEigenVectors( bcSide, *M_solution_tn, eigenvalues, leftEigenvector1, leftEigenvector2 );

            switch ( bcSide )
            {
            case OneD_left:
                switch ( bcOutputType )
                {
                case OneD_Q: // dQ_L/dP_L
                    JacobianCoefficient = leftEigenvector2[0] / leftEigenvector2[1]
                                          * M_physics->dAdP( M_solver->boundaryValue( *M_solution, OneD_P, OneD_left ), 0 );
                    break;
                case OneD_P: // dP_L/dQ_L
                    JacobianCoefficient = leftEigenvector2[1] / leftEigenvector2[0]
                                          * M_physics->dPdA( M_solver->boundaryValue( *M_solution, OneD_A, OneD_left ), 0 );
                    break;
                default:
                    std::cout << "Warning: bcType \"" << bcOutputType << "\"not available!" << std::endl;
                }
                break;
            case OneD_right:
                switch ( bcOutputType )
                {
                case OneD_Q: // dQ_R/dP_R
                    JacobianCoefficient = -leftEigenvector1[0] / leftEigenvector1[1]
                                          * M_physics->dAdP( M_solver->boundaryValue( *M_solution, OneD_P, OneD_right ), M_data->NumberOfElements() );
                    break;
                case OneD_P: // dP_R/dQ_R
                    JacobianCoefficient = -leftEigenvector1[1] / leftEigenvector1[0]
                                          * M_physics->dPdA( M_solver->boundaryValue( *M_solution, OneD_A, OneD_right ), M_data->NumberOfElements() );
                    break;
                default:
                    std::cout << "Warning: bcType \"" << bcOutputType << "\"not available!" << std::endl;
                }
                break;
            default:
                std::cout << "Warning: bcSide \"" << bcSide << "\" not available!" << std::endl;
            }

            break;
        }

    return jacobianCoefficient;
}

#endif

} // Namespace LifeV
