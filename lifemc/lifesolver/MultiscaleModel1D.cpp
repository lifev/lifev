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
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModel1D::MultiscaleModel1D() :
        multiscaleModel_Type           (),
#ifdef HAVE_HDF5
        M_exporter                     ( new IOFile_Type() ),
        M_importer                     ( new IOFile_Type() ),
        M_exporterMesh                 ( new mesh_Type() ),
#endif
#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
        M_linearBC                     ( new bc_Type() ),
        M_linearSolution               ( new solution_Type() ),
        M_bcPreviousTimeSteps          (),
        M_bcBaseDelta                  (),
        M_bcDelta                      (),
        M_bcDeltaType                  (),
        M_bcDeltaSide                  (),
#endif
        M_data                         ( new data_Type() ),
        M_bc                           ( new bcInterface_Type() ),
        M_physics                      (),
        M_flux                         (),
        M_source                       (),
        M_solver                       ( new solver_Type() ),
        M_linearSolver                 (),
        M_FESpace                      (),
        M_solution_tn                  ( new solution_Type() ),
        M_solution                     ( new solution_Type() )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::MultiscaleModel1D() \n";
#endif

    M_type = OneDimensional;

    //Define the maps of the OneDimensionalModel objects
    OneDimensional::mapsDefinition();

    //Register the objects
//!\todo pass a std::string to the factories
//     factoryOneDimensionalPhysics_Type::instance().registerProduct( "OneD_LinearPhysics",    &createOneDimensionalPhysicsLinear );
//     factoryOneDimensionalPhysics_Type::instance().registerProduct( "OneD_NonLinearPhysics", &createOneDimensionalPhysicsNonLinear );

//     factoryOneDimensionalFlux_Type::instance().registerProduct(    "OneD_LinearFlux",       &createOneDimensionalFluxLinear );
//     factoryOneDimensionalFlux_Type::instance().registerProduct(    "OneD_NonLinearFlux",    &createOneDimensionalFluxNonLinear );

//     factoryOneDimensionalSource_Type::instance().registerProduct(  "OneD_LinearSource",     &createOneDimensionalSourceLinear );
//     factoryOneDimensionalSource_Type::instance().registerProduct(  "OneD_NonLinearSource",  &createOneDimensionalSourceNonLinear );

    physics_Type::factoryPhysics_Type::instance().registerProduct( OneDimensional::LinearPhysics,    &createOneDimensionalPhysicsLinear );
    physics_Type::factoryPhysics_Type::instance().registerProduct( OneDimensional::NonLinearPhysics, &createOneDimensionalPhysicsNonLinear );

    flux_Type::factoryFlux_Type::instance().registerProduct(       OneDimensional::LinearFlux,       &createOneDimensionalFluxLinear );
    flux_Type::factoryFlux_Type::instance().registerProduct(       OneDimensional::NonLinearFlux,    &createOneDimensionalFluxNonLinear );

    source_Type::factorySource_Type::instance().registerProduct(   OneDimensional::LinearSource,     &createOneDimensionalSourceLinear );
    source_Type::factorySource_Type::instance().registerProduct(   OneDimensional::NonLinearSource,  &createOneDimensionalSourceNonLinear );
}

// ===================================================
// MultiScale PhysicalModel Virtual Methods
// ===================================================
void
MultiscaleModel1D::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::SetupData( ) \n";
#endif

    // Preliminary setup of the communicator
#ifdef EPETRA_MPI
    MPI_Comm LocalComm;
    MPI_Comm_split( ( dynamic_cast<Epetra_MpiComm*> ( &(*M_comm) ) )->Comm(), M_comm->MyPID(), M_comm->MyPID(), &LocalComm );
    M_comm.reset( new Epetra_MpiComm( LocalComm ) );
#else
    M_comm.reset( new Epetra_SerialComm() );
#endif

    multiscaleModel_Type::setupData( fileName );

    GetPot dataFile( fileName );

    M_data->setup( dataFile );
    if ( M_globalData.get() )
        setupGlobalData( fileName );

    //1D Model Physics
    //!\todo pass a std::string to the factories
    //M_physics = physicsPtr_Type( factoryOneDimensionalPhysics_Type::instance().createObject( "M_data->physicsType()" ) );
    M_physics = physicsPtr_Type( physics_Type::factoryPhysics_Type::instance().createObject( M_data->physicsType() ) );
    M_physics->setData( M_data );

    //1D Model Flux
    //!\todo pass a std::string to the factories
    //M_flux = fluxPtr_Type( factoryOneDimensionalFlux_Type::instance().createObject( "M_data->fluxType()" ) );
    M_flux = fluxPtr_Type( flux_Type::factoryFlux_Type::instance().createObject( M_data->fluxType() ) );
    M_flux->setPhysics( M_physics );

    //1D Model Source
    //!\todo pass a std::string to the factories
    //M_source = sourcePtr_Type( factoryOneDimensionalSource_Type::instance().createObject( "M_data->sourceType()" ) );
    M_source = sourcePtr_Type( source_Type::factorySource_Type::instance().createObject( M_data->sourceType() ) );
    M_source->setPhysics( M_physics );

    //Linear Solver
    M_linearSolver.reset( new linearSolver_Type( M_comm ) );
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
    M_data->setPostprocessingDirectory( multiscaleProblemFolder );
    M_data->setPostprocessingFile( "Step_" + number2string( multiscaleProblemStep ) + "_Model_" + number2string( M_ID ) );

#ifdef HAVE_HDF5
    M_exporter->setDataFromGetPot( dataFile );
    M_exporter->setPrefix( "Step_" + number2string( multiscaleProblemStep ) + "_Model_" + number2string( M_ID ) );
    M_exporter->setPostDir( multiscaleProblemFolder );

    M_exporterMesh->setup( M_data->length(), M_data->numberOfElements() );

    M_importer->setDataFromGetPot( dataFile );
    M_importer->setPrefix( "Step_" + number2string( multiscaleProblemStep - 1 ) + "_Model_" + number2string( M_ID ) );
    M_importer->setPostDir( multiscaleProblemFolder );
#endif

}

void
MultiscaleModel1D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::SetupProblem() \n";
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
MultiscaleModel1D::buildSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::BuildSystem() \n";
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
MultiscaleModel1D::updateSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::UpdateSystem() \n";
#endif

    // Update previous solution
    updateSolution( *M_solution, *M_solution_tn );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
        updateLinearModel();
#endif

}

void
MultiscaleModel1D::solveSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::SolveSystem() \n";
#endif

    solve( *M_bc->handler(), *M_solution );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
        updateLinearBC( *M_solution );
#endif

}

void
MultiscaleModel1D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::SaveSolution() \n";
#endif

#ifdef HAVE_HDF5
    M_exporter->postProcess( M_data->dataTime()->time() );

    if ( M_data->dataTime()->isLastTimeStep() )
        M_exporter->closeFile();
#endif

#ifdef HAVE_MATLAB_POSTPROCESSING
    //Matlab post-processing
    M_solver->postProcess( *M_solution );
#endif

}

void
MultiscaleModel1D::showMe()
{
    if ( M_displayer->isLeader() )
    {
        multiscaleModel_Type::showMe();

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
MultiscaleModel1D::setupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::SetupLinearModel( ) \n";
#endif

    // Define bcFunction for linear problem
    M_bcBaseDelta.setFunction( boost::bind( &MultiscaleModel1D::bcFunctionDelta, this, _1 ) );

    // The linear BCHandler is a copy of the original BCHandler with the LinearSolution instead of the true solution
    //M_LinearBC.reset( new bc_Type( *M_bc->handler() ) ); // COPY CONSTRUCTOR NOT WORKING

    //Set left and right BC + default BC
    M_linearBC->setBC( OneDimensional::left, OneDimensional::first, M_bc->handler()->bc( OneDimensional::left )->type( OneDimensional::first ),
                       M_bc->handler()->bc( OneDimensional::left )->bcFunction( OneDimensional::first ) );

    M_linearBC->setBC( OneDimensional::right, OneDimensional::first, M_bc->handler()->bc( OneDimensional::right )->type( OneDimensional::first ),
                       M_bc->handler()->bc( OneDimensional::right )->bcFunction( OneDimensional::first ) );

    M_linearBC->setDefaultBC();

    // Solution for the linear problem (this does not change anything in the solver)
    M_solver->setupSolution( *M_linearSolution );
    M_linearBC->setSolution( M_linearSolution );
    M_linearBC->setFluxSource( M_flux, M_source );
}

void
MultiscaleModel1D::updateLinearModel()
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
MultiscaleModel1D::solveLinearModel( bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::SolveLinearModel() \n";
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
Real
MultiscaleModel1D::boundaryStress( const bcFlag_Type& flag, const stress_Type& stressType ) const
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
        std::cout << "ERROR: Invalid stress type [" << Enum2String( stressType, multiscaleStressesMap ) << "]" << std::endl;

        return 0.0;
    }
    }
}

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE

Real
MultiscaleModel1D::boundaryDeltaFlowRate( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    bcSide_Type bcSide = flagConverter( flag );

    solveLinearModel( solveLinearSystem );

    Real Q      = M_solver->boundaryValue( *M_solution, OneDimensional::Q, bcSide );
    Real Qdelta = M_solver->boundaryValue( *M_linearSolution, OneDimensional::Q, bcSide );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::GetBoundaryDeltaFlowRate( flag, solveLinearSystem ) \n";
    Debug( 8130 ) << "Q:          " << Q << "\n";
    Debug( 8130 ) << "Qdelta:     " << Qdelta << "\n";
#endif

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA

    if ( M_bcDeltaType == OneDimensional::A )
    {
        // dQ/dP
        return ( (Qdelta - Q) / M_bcDelta ) * M_physics->dAdP( M_solver->boundaryValue( *M_solution, OneDimensional::P, bcSide ), 0 );
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
MultiscaleModel1D::boundaryDeltaPressure( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    bcSide_Type bcSide = flagConverter( flag );

    solveLinearModel( solveLinearSystem );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA

    Real A      = M_solver->boundaryValue( *M_solution, OneDimensional::A, bcSide );
    Real Adelta = M_solver->boundaryValue( *M_linearSolution, OneDimensional::A, bcSide );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::GetBoundaryDeltaPressure( flag, solveLinearSystem ) \n";
    Debug( 8130 ) << "A:          " << A <<  "\n";
    Debug( 8130 ) << "Adelta:     " << Adelta <<  "\n";
#endif

    if ( M_bcDeltaType == OneDimensional::A )
    {
        // dP/dP
        return ( (Adelta - A) / M_bcDelta ) * M_physics->dPdA( M_solver->boundaryValue( *M_solution, OneDimensional::A, bcSide ), 0 )
               * M_physics->dAdP( M_solver->boundaryValue( *M_solution, OneDimensional::P, bcSide ), 0 );
    }
    else
    {
        // dP/dQ
        return ( (Adelta - A) / M_bcDelta ) * M_physics->dPdA( M_solver->boundaryValue( *M_solution, OneDimensional::A, bcSide ), 0 );
    }

#else

    Real P      = M_solver->boundaryValue( *M_solution, OneDimensional::P, bcSide );
    Real Pdelta = M_solver->boundaryValue( *M_linearSolution, OneDimensional::P, bcSide );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::GetBoundaryDeltaPressure( flag, solveLinearSystem ) \n";
    Debug( 8130 ) << "P:          " << P <<  "\n";
    Debug( 8130 ) << "Pdelta:     " << Pdelta <<  "\n";
#endif

    return (Pdelta - P) / M_bcDelta;

#endif

}

#else

Real
MultiscaleModel1D::boundaryDeltaFlowRate( const bcFlag_Type& flag, bool& /*solveLinearSystem*/ )
{
    return tangentProblem( flagConverter( flag ), OneDimensional::Q );
}

Real
MultiscaleModel1D::boundaryDeltaPressure( const bcFlag_Type& flag, bool& /*solveLinearSystem*/ )
{
    return tangentProblem( flagConverter( flag ), OneDimensional::P );
}

#endif

Real
MultiscaleModel1D::boundaryDeltaDynamicPressure( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    // Untested
    return boundaryDensity( flag ) * boundaryDeltaFlowRate( flag, solveLinearSystem ) * boundaryFlowRate( flag ) / ( boundaryArea( flag ) * boundaryArea( flag ) );
}

Real
MultiscaleModel1D::boundaryDeltaStress( const bcFlag_Type& flag, bool& solveLinearSystem, const stress_Type& stressType )
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
        std::cout << "ERROR: Invalid stress type [" << Enum2String( stressType, multiscaleStressesMap ) << "]" << std::endl;

        return 0.0;
    }
    }
}

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleModel1D::setupGlobalData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::SetupGlobalData( fileName ) \n";
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
MultiscaleModel1D::setupFESpace()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::SetupFEspace() \n";
#endif

    //Transform mesh
    boost::array< Real, NDIM > NullTransformation;
    NullTransformation[0] = 0.;
    NullTransformation[1] = 0.;
    NullTransformation[2] = 0.;

    //The real mesh can be only scaled due to OneDimensionalSolver conventions
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

    M_FESpace.reset( new feSpace_Type( M_data->mesh(), *refFE, *qR, *bdQr, 1, M_comm ) );
    M_solver->setFESpace( M_FESpace );
}

void
MultiscaleModel1D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::InitializeSolution() \n";
#endif

    if ( multiscaleProblemStep > 0 )
    {
        M_importer->setMeshProcId( M_exporterMesh, M_comm->MyPID() );

        //M_exporter->addVariable( ExporterData::Scalar, "Solid Area",      (*M_Solution)["A"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        M_importer->addVariable( ExporterData::Scalar, "Area ratio",      (*M_solution)["A/A0-1"], static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        M_importer->addVariable( ExporterData::Scalar, "Fluid Flow Rate", (*M_solution)["Q"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        //M_importer->addVariable( ExporterData::Scalar, "W1",              (*M_Solution)["W1"],   static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        //M_importer->addVariable( ExporterData::Scalar, "W2",              (*M_Solution)["W2"],   static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        M_importer->addVariable( ExporterData::Scalar, "Fluid Pressure",  (*M_solution)["P"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );

        // Import
        M_exporter->setStartIndex( M_importer->importFromTime( M_data->dataTime()->initialTime() ) + 1 );

        // Compute A from AreaRatio
        M_solver->computeArea( *M_solution );

        // Compute W1 and W2 from A and Q
        M_solver->computeW1W2( *M_solution );
    }
    else
        M_solver->initialize( *M_solution );
}

void
MultiscaleModel1D::updateSolution( const solution_Type& solution1, solution_Type& solution2 )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::UpdateSolution( solution1, solution2 ) \n";
#endif

    // Here we make a true copy (not a copy of the shared_ptr, but a copy of its content)
    for ( solutionConstIterator_Type i = solution1.begin() ; i != solution1.end() ; ++i )
        *solution2[i->first] = *i->second;
}

void
MultiscaleModel1D::solve( bc_Type& bc, solution_Type& solution, const std::string& solverType )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::Solve() \n";
#endif

    // Re-initialize solution
    updateSolution( *M_solution_tn, solution );

    // Subiterate to respect CFL
    UInt SubiterationNumber(1);
    Real timeStep = M_data->dataTime()->timeStep();

    Real CFL = M_solver->computeCFL( solution, M_data->dataTime()->timeStep() );
    if ( CFL > M_data->CFLmax() )
    {
        SubiterationNumber = std::ceil( CFL / M_data->CFLmax() );
        timeStep /= SubiterationNumber;
    }

    if ( M_displayer->isLeader() )
        std::cout << solverType << "  CFL                                      " << CFL*timeStep/M_data->dataTime()->timeStep() << std::endl;

    for ( UInt i(1) ; i <= SubiterationNumber ; ++i )
    {
        if ( M_displayer->isLeader() )
        {
            std::cout << solverType << "  Subiteration                             " << i << "/" << SubiterationNumber << std::endl;
            std::cout << solverType << "  Time                                     " <<  M_data->dataTime()->previousTime() + i*timeStep << std::endl;
        }
        //bc.updateOperatorVariables();

        M_solver->updateRHS( solution, timeStep );
        M_solver->iterate( bc, solution, M_data->dataTime()->previousTime() + i*timeStep, timeStep );
    }
}

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE

void
MultiscaleModel1D::createLinearBC()
{
    // Allocating the correct space
    M_bcPreviousTimeSteps.reserve( std::max( M_couplings[0]->timeInterpolationOrder(), M_couplings[1]->timeInterpolationOrder() ) );

    // Create bcSide map
    std::map< bcSide_Type, std::map< bcType_Type, Real > > bcSideMap;
    M_bcPreviousTimeSteps.push_back( bcSideMap );

    // Create bcType map
    std::map< bcType_Type, Real > bcTypeMap;
    M_bcPreviousTimeSteps[0][OneDimensional::left]  = bcTypeMap;
    M_bcPreviousTimeSteps[0][OneDimensional::right] = bcTypeMap;
}

void
MultiscaleModel1D::updateLinearBC( const solution_Type& solution )
{
    M_bcPreviousTimeSteps[0][OneDimensional::left][OneDimensional::A]  = M_solver->boundaryValue( solution, OneDimensional::A, OneDimensional::left );
    M_bcPreviousTimeSteps[0][OneDimensional::left][OneDimensional::P]  = M_solver->boundaryValue( solution, OneDimensional::P, OneDimensional::left );
    M_bcPreviousTimeSteps[0][OneDimensional::left][OneDimensional::Q]  = M_solver->boundaryValue( solution, OneDimensional::Q, OneDimensional::left );
    M_bcPreviousTimeSteps[0][OneDimensional::right][OneDimensional::A] = M_solver->boundaryValue( solution, OneDimensional::A, OneDimensional::right );
    M_bcPreviousTimeSteps[0][OneDimensional::right][OneDimensional::P] = M_solver->boundaryValue( solution, OneDimensional::P, OneDimensional::right );
    M_bcPreviousTimeSteps[0][OneDimensional::right][OneDimensional::Q] = M_solver->boundaryValue( solution, OneDimensional::Q, OneDimensional::right );
}

void
MultiscaleModel1D::imposePerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::Perturbation() \n";
#endif

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            // Find the side to perturb and apply the perturbation
            M_bcDeltaSide = flagConverter( ( *i )->flag( ( *i )->modelGlobalToLocalID( M_ID ) ) );
            M_linearBC->bc( M_bcDeltaSide )->setBCFunction( OneDimensional::first, M_bcBaseDelta );

            // Compute the range
            M_bcDeltaType = M_linearBC->bc( M_bcDeltaSide )->type( OneDimensional::first );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA
            // We replace pressure BC with area BC for the perturbed problem
            if ( M_bcDeltaType == OneDimensional::P )
            {
                M_linearBC->bc( M_bcDeltaSide )->setType( OneDimensional::first, OneDimensional::A );
                M_bcDeltaType = OneDimensional::A;
            }
#endif

            //M_BCDelta = ( *i )->residual()[( *i )->perturbedCoupling()] * 10000;
            //M_BCDelta = ( M_BCDelta[1] - M_BCDelta[0] ) / 100;

            //if ( std::abs( M_BCDelta ) < 1e-6 || std::abs( M_BCDelta ) > 1e6 )
            switch ( M_bcDeltaType )
            {
            case OneDimensional::A:

                M_bcDelta = M_data->jacobianPerturbationArea();
                break;

            case OneDimensional::Q:

                M_bcDelta = M_data->jacobianPerturbationFlowRate();

                break;

            case OneDimensional::P:

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
MultiscaleModel1D::resetPerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::ResetPerturbation() \n";
#endif

    M_linearBC->bc( M_bcDeltaSide )->setBCFunction( OneDimensional::first, M_bc->handler()->bc( M_bcDeltaSide )->bcFunction( OneDimensional::first ) );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA
    // Restoring the original BC
    if ( M_bcDeltaType == OneDimensional::A )
        M_linearBC->bc( M_bcDeltaSide )->setType( OneDimensional::first, OneDimensional::P );
#endif

}

Real
MultiscaleModel1D::bcFunctionDelta( const Real& t )
{
    // Lagrange interpolation
    Real bcValue(0);
    Real base(1);

    // Time container for interpolation
    std::vector< Real > timeContainer( M_bcPreviousTimeSteps.size(), 0 );
    for ( UInt i(0) ; i < M_bcPreviousTimeSteps.size() ; ++i )
        timeContainer[i] = M_globalData->dataTime()->time() - i * M_globalData->dataTime()->timeStep();

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
MultiscaleModel1D::tangentProblem( const bcSide_Type& bcOutputSide, const bcType_Type& bcOutputType )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::TangentProblem( bcOutputSide, bcOutputType ) \n";
#endif

    Real jacobianCoefficient(0);

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            // Find the perturbed side
            bcSide_Type bcSide = flagConverter( ( *i )->flag( ( *i )->modelGlobalToLocalID( M_ID ) ) );

            // Perturbation has no effect on the other sides (which also means that dQ/dQ and dP/dP are always zero)
            if ( bcSide != bcOutputSide )
                break;

            // Compute the eigenvectors
            Container2D_Type eigenvalues, leftEigenvector1, leftEigenvector2;
            M_solver->BoundaryEigenValuesEigenVectors( bcSide, *M_solution_tn, eigenvalues, leftEigenvector1, leftEigenvector2 );

            switch ( bcSide )
            {
            case OneDimensional::left:
                switch ( bcOutputType )
                {
                case OneDimensional::Q: // dQ_L/dP_L
                    JacobianCoefficient = leftEigenvector2[0] / leftEigenvector2[1]
                                          * M_physics->dAdP( M_solver->boundaryValue( *M_solution, OneDimensional::P, OneDimensional::left ), 0 );
                    break;
                case OneDimensional::P: // dP_L/dQ_L
                    JacobianCoefficient = leftEigenvector2[1] / leftEigenvector2[0]
                                          * M_physics->dPdA( M_solver->boundaryValue( *M_solution, OneDimensional::A, OneDimensional::left ), 0 );
                    break;
                default:
                    std::cout << "Warning: bcType \"" << bcOutputType << "\"not available!" << std::endl;
                }
                break;
            case OneDimensional::right:
                switch ( bcOutputType )
                {
                case OneDimensional::Q: // dQ_R/dP_R
                    JacobianCoefficient = -leftEigenvector1[0] / leftEigenvector1[1]
                                          * M_physics->dAdP( M_solver->boundaryValue( *M_solution, OneDimensional::P, OneDimensional::right ), M_data->NumberOfElements() );
                    break;
                case OneDimensional::P: // dP_R/dQ_R
                    JacobianCoefficient = -leftEigenvector1[1] / leftEigenvector1[0]
                                          * M_physics->dPdA( M_solver->boundaryValue( *M_solution, OneDimensional::A, OneDimensional::right ), M_data->NumberOfElements() );
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

} // Namespace multiscale
} // Namespace LifeV
