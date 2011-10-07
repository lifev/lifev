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
 *  @brief File containing the Multiscale Model 1D
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
        MultiscaleInterfaceFluid       (),
#ifdef HAVE_HDF5
        M_exporter                     ( new IOFile_Type() ),
        M_importer                     ( new IOFile_Type() ),
        M_exporterMesh                 ( new mesh_Type() ),
        M_exporterSolution             ( new solution_Type() ),
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
        M_linearViscoelasticSolver     (),
        M_feSpace                      (),
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
    physics_Type::factoryPhysics_Type::instance().registerProduct( OneDimensional::LinearPhysics,    &createOneDimensionalPhysicsLinear );
    physics_Type::factoryPhysics_Type::instance().registerProduct( OneDimensional::NonLinearPhysics, &createOneDimensionalPhysicsNonLinear );

    flux_Type::factoryFlux_Type::instance().registerProduct(       OneDimensional::LinearFlux,       &createOneDimensionalFluxLinear );
    flux_Type::factoryFlux_Type::instance().registerProduct(       OneDimensional::NonLinearFlux,    &createOneDimensionalFluxNonLinear );

    source_Type::factorySource_Type::instance().registerProduct(   OneDimensional::LinearSource,     &createOneDimensionalSourceLinear );
    source_Type::factorySource_Type::instance().registerProduct(   OneDimensional::NonLinearSource,  &createOneDimensionalSourceNonLinear );
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModel1D::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::setupData( fileName ) \n";
#endif

    multiscaleModel_Type::setupData( fileName );

    GetPot dataFile( fileName );

    M_data->setup( dataFile );
    if ( M_globalData.get() )
        setupGlobalData( fileName );

    //1D Model Physics
    M_physics = physicsPtr_Type( physics_Type::factoryPhysics_Type::instance().createObject( M_data->physicsType(), OneDimensional::physicsMap ) );
    M_physics->setData( M_data );

    //1D Model Flux
    M_flux = fluxPtr_Type( flux_Type::factoryFlux_Type::instance().createObject( M_data->fluxType(), OneDimensional::fluxMap ) );
    M_flux->setPhysics( M_physics );

    //1D Model Source
    M_source = sourcePtr_Type( source_Type::factorySource_Type::instance().createObject( M_data->sourceType(), OneDimensional::sourceMap ) );
    M_source->setPhysics( M_physics );

    //Linear Solver
    M_linearSolver.reset( new linearSolver_Type( M_comm ) );
    //M_linearSolver->setupPreconditioner( dataFile, "1D_Model/prec" );
    M_linearSolver->setDataFromGetPot( dataFile, "1D_Model/solver" );
    M_linearSolver->setParameter( "Verbose", false );
    M_linearSolver->setParameters();

    //Linear Viscoelastic Solver
    if ( M_data->viscoelasticWall() )
    {
        M_linearViscoelasticSolver.reset( new linearSolver_Type( M_comm ) );
        M_linearViscoelasticSolver->setParametersList( M_linearSolver->parametersList() );
        M_linearViscoelasticSolver->setParameters();
    }

    //1D Model Solver
    M_solver->setCommunicator( M_comm );
    M_solver->setProblem( M_physics, M_flux, M_source );
    M_solver->setLinearSolver( M_linearSolver );
    M_solver->setLinearViscoelasticSolver( M_linearViscoelasticSolver );

    //BC - We need to create the BCHandler before using it
    M_bc->createHandler();
    //M_bc->fillHandler( fileName, "1D_Model" );

    //Exporters
    M_data->setPostprocessingDirectory( multiscaleProblemFolder );
    M_data->setPostprocessingFile( multiscaleProblemPrefix + "_Model_" + number2string( M_ID ) + "_" + number2string( multiscaleProblemStep ) );

#ifdef HAVE_HDF5
    M_exporter->setDataFromGetPot( dataFile );
    M_exporter->setPrefix( multiscaleProblemPrefix + "_Model_" + number2string( M_ID ) + "_" + number2string( multiscaleProblemStep ) );
    M_exporter->setPostDir( multiscaleProblemFolder );
    M_exporterMesh->setup( M_data->length(), M_data->numberOfElements() );

    M_importer->setDataFromGetPot( dataFile );
    M_importer->setPrefix( multiscaleProblemPrefix + "_Model_" + number2string( M_ID ) + "_" + number2string( multiscaleProblemStep -1 ) );
    M_importer->setPostDir( multiscaleProblemFolder );
#endif

}

void
MultiscaleModel1D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::setupModel() \n";
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

    //Post-processing
#ifdef HAVE_HDF5
    M_exporter->setMeshProcId( M_exporterMesh, M_comm->MyPID() );

    MapEpetra map( M_feSpace->refFE(), *M_exporterMesh, M_comm );
    M_solver->setupSolution( *M_exporterSolution, map, true );

    M_exporter->addVariable( IOData_Type::ScalarField, "Area ratio (fluid)", M_feSpace, (*M_exporterSolution)["AoverA0minus1"], static_cast <UInt> ( 0 ) );
    M_exporter->addVariable( IOData_Type::ScalarField, "Flow rate (fluid)",  M_feSpace, (*M_exporterSolution)["Q"],    static_cast <UInt> ( 0 ) );
    //M_exporter->addVariable( IOData_Type::ScalarField, "W1",               M_feSpace, (*M_exporterSolution)["W1"],   static_cast <UInt> ( 0 ), M_feSpace->dof().numTotalDof() );
    //M_exporter->addVariable( IOData_Type::ScalarField, "W2",               M_feSpace, (*M_exporterSolution)["W2"],   static_cast <UInt> ( 0 ), M_feSpace->dof().numTotalDof() );
    M_exporter->addVariable( IOData_Type::ScalarField, "Pressure (fluid)",   M_feSpace, (*M_exporterSolution)["P"],    static_cast <UInt> ( 0 ) );
#endif

#ifdef HAVE_MATLAB_POSTPROCESSING
    M_solver->resetOutput( *M_exporterSolution );
#endif

    //Setup solution
    initializeSolution();

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
    {
        createLinearBC();
        updateLinearBC( *M_solution );
        setupLinearModel();

        // Initialize the linear solution
        copySolution( *M_solution, *M_linearSolution );
    }
#endif

}

void
MultiscaleModel1D::buildModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::buildModel() \n";
#endif

    // Display data
//    if ( M_comm->MyPID() == 0 )
//        M_data->showMe();

    M_solver->buildConstantMatrices();

    // Update previous solution
    copySolution( *M_solution, *M_solution_tn );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
        updateLinearModel();
#endif

}

void
MultiscaleModel1D::updateModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::updateModel() \n";
#endif

    // Update previous solution
    copySolution( *M_solution, *M_solution_tn );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
        updateLinearModel();
#endif

}

void
MultiscaleModel1D::solveModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::solveModel() \n";
#endif

    displayModelStatus( "Solve" );
    solve( *M_bc->handler(), *M_solution );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
        updateLinearBC( *M_solution );
#endif

}

void
MultiscaleModel1D::updateSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::updateSolution() \n";
#endif

    // Update exporter solution removing ghost nodes
    copySolution( *M_solution, *M_exporterSolution );
}

void
MultiscaleModel1D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::saveSolution() \n";
#endif

#ifdef HAVE_HDF5
    M_exporter->postProcess( M_data->dataTime()->time() );

    if ( M_data->dataTime()->isLastTimeStep() )
        M_exporter->closeFile();
#endif

#ifdef HAVE_MATLAB_POSTPROCESSING
    //Matlab post-processing
    M_solver->postProcess( *M_exporterSolution, M_data->dataTime()->time() );
#endif

}

void
MultiscaleModel1D::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        multiscaleModel_Type::showMe();

        std::cout << "FE order            = " << "P1" << std::endl
                  << "DOF                 = " << M_data->mesh()->numPoints() << std::endl << std::endl;

        std::cout << "maxH                = " << M_data->mesh()->maxH() << std::endl
                  << "meanH               = " << M_data->mesh()->meanH() << std::endl << std::endl;
    }
}

Real
MultiscaleModel1D::checkSolution() const
{
    return (*M_solution)["AoverA0minus1"]->norm2() + (*M_solution)["Q"]->norm2() + (*M_solution)["P"]->norm2();
}

// ===================================================
// MultiscaleInterfaceFluid Methods
// ===================================================
void
MultiscaleModel1D::imposeBoundaryFlowRate( const bcFlag_Type& flag, const function_Type& function )
{
    OneDimensionalFunction base;
    base.setFunction( boost::bind( function, _1, _1, _1, _1, _1 ) );

    M_bc->handler()->setBC( flagConverter( flag ), OneDimensional::first, OneDimensional::Q, base );
}

void
MultiscaleModel1D::imposeBoundaryStress( const bcFlag_Type& flag, const function_Type& function )
{
    OneDimensionalFunction base;
    base.setFunction( boost::bind( function, _1, _1, _1, _1, _1 ) );

    M_bc->handler()->setBC( flagConverter( flag ), OneDimensional::first, OneDimensional::S, base );
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
    Debug( 8130 ) << "MultiscaleModel1D::boundaryDeltaFlowRate( flag, solveLinearSystem ) \n";
    Debug( 8130 ) << "Q:          " << Q << "\n";
    Debug( 8130 ) << "Qdelta:     " << Qdelta << "\n";
#endif

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA

    if ( M_bcDeltaType == OneDimensional::A )
    {
        // dQ/dP
        return ( (Qdelta - Q) / M_bcDelta ) * M_physics->dAdP( M_solver->boundaryValue( *M_solution, OneDimensional::P, bcSide ), M_data->dataTime()->timeStep(), 0 );
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
MultiscaleModel1D::boundaryDeltaStress( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    bcSide_Type bcSide = flagConverter( flag );

    solveLinearModel( solveLinearSystem );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA

    Real A      = M_solver->boundaryValue( *M_solution, OneDimensional::A, bcSide );
    Real Adelta = M_solver->boundaryValue( *M_linearSolution, OneDimensional::A, bcSide );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::boundaryDeltaStress( flag, solveLinearSystem ) \n";
    Debug( 8130 ) << "A:          " << A <<  "\n";
    Debug( 8130 ) << "Adelta:     " << Adelta <<  "\n";
#endif

    if ( M_bcDeltaType == OneDimensional::A )
    {
        // dP/dP
        return ( (Adelta - A) / M_bcDelta ) * M_physics->dPdA( M_solver->boundaryValue( *M_solution, OneDimensional::A, bcSide ), M_data->dataTime()->timeStep(), 0 )
               * M_physics->dAdP( M_solver->boundaryValue( *M_solution, OneDimensional::P, bcSide ), M_data->dataTime()->timeStep(), 0 );
    }
    else
    {
        // dP/dQ
        return ( (Adelta - A) / M_bcDelta ) * M_physics->dPdA( M_solver->boundaryValue( *M_solution, OneDimensional::A, bcSide ), M_data->dataTime()->timeStep(), 0 );
    }

#else

    Real S      = M_solver->boundaryValue( *M_solution, OneDimensional::S, bcSide );
    Real Sdelta = M_solver->boundaryValue( *M_linearSolution, OneDimensional::S, bcSide );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::boundaryDeltaStress( flag, solveLinearSystem ) \n";
    Debug( 8130 ) << "S:          " << S <<  "\n";
    Debug( 8130 ) << "Sdelta:     " << Sdelta <<  "\n";
#endif

    return (Sdelta - S) / M_bcDelta;
#endif
}

#else

Real
MultiscaleModel1D::boundaryDeltaFlowRate( const bcFlag_Type& flag, bool& /*solveLinearSystem*/ )
{
    return tangentProblem( flagConverter( flag ), OneDimensional::Q );
}

Real
MultiscaleModel1D::boundaryDeltaStress( const bcFlag_Type& flag, bool& /*solveLinearSystem*/ )
{
    return tangentProblem( flagConverter( flag ), OneDimensional::S );
}

#endif

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleModel1D::setupGlobalData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::setupGlobalData( fileName ) \n";
#endif

    GetPot dataFile( fileName );

    //Global data time
    M_data->setTimeData( M_globalData->dataTime() );

    //Global physical quantities
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/density" ) )
        M_data->setDensity( M_globalData->fluidDensity() );
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/viscosity" ) )
        M_data->setViscosity( M_globalData->fluidViscosity() );
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/venousPressure" ) )
        M_data->setVenousPressure( M_globalData->fluidVenousPressure() );
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/externalPressure" ) )
        M_data->setExternalPressure( M_globalData->solidExternalPressure() );
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/densityWall" ) )
        M_data->setDensityWall( M_globalData->solidDensity() );
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/poisson" ) )
        M_data->setPoisson( M_globalData->solidPoissonCoefficient() );
    if ( !dataFile.checkVariable( "1D_Model/PhysicalParameters/young" ) )
        M_data->setYoung( M_globalData->solidYoungModulus() );

    //After changing some parameters we need to update the coefficients
    M_data->updateCoefficients();
}

void
MultiscaleModel1D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::initializeSolution() \n";
#endif

    if ( multiscaleProblemStep > 0 )
    {
#ifdef HAVE_HDF5
        M_importer->setMeshProcId( M_exporterMesh, M_comm->MyPID() );

        M_importer->addVariable( IOData_Type::ScalarField, "Area ratio (fluid)", M_feSpace, (*M_exporterSolution)["AoverA0minus1"], static_cast <UInt> ( 0 ) );
        M_importer->addVariable( IOData_Type::ScalarField, "Flow rate (fluid)",  M_feSpace, (*M_exporterSolution)["Q"],      static_cast <UInt> ( 0 ) );
//        M_importer->addVariable( IOData_Type::ScalarField, "W1",               M_feSpace, (*M_exporterSolution)["W1"],     static_cast <UInt> ( 0 ), M_feSpace->dof().numTotalDof() );
//        M_importer->addVariable( IOData_Type::ScalarField, "W2",               M_feSpace, (*M_exporterSolution)["W2"],     static_cast <UInt> ( 0 ), M_feSpace->dof().numTotalDof() );
        M_importer->addVariable( IOData_Type::ScalarField, "Pressure (fluid)",   M_feSpace, (*M_exporterSolution)["P"],      static_cast <UInt> ( 0 ) );

        // Import
        M_exporter->setTimeIndex( M_importer->importFromTime( M_data->dataTime()->initialTime() ) + 1 );

        M_importer->closeFile();
#else
        std::cout << "!!! ERROR: Importer not implemented for this filter !!!" << std::endl;
#endif

        // Copy the imported solution to the problem solution container
        copySolution( *M_exporterSolution, *M_solution );

        // Compute A from AreaRatio
        M_solver->computeArea( *M_solution );

        // Compute W1 and W2 from A and Q
        M_solver->computeW1W2( *M_solution );
    }
    else
    {
        M_solver->initialize( *M_solution );
        copySolution( *M_solution, *M_exporterSolution );
    }
}

void
MultiscaleModel1D::setupFESpace()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::setupFEspace() \n";
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
    const ReferenceFE*    refFE = &feSegP1;
    const QuadratureRule* qR    = &quadRuleSeg3pt;
    const QuadratureRule* bdQr  = &quadRuleSeg1pt;

//    const RefFE*    refFE = &feSegP2;
//    const QuadRule* qR    = &quadRuleSeg3pt;
//    const QuadRule* bdQr  = &quadRuleSeg1pt;

    M_feSpace.reset( new feSpace_Type( M_data->mesh(), *refFE, *qR, *bdQr, 1, M_comm ) );
    M_solver->setFESpace( M_feSpace );
}

void
MultiscaleModel1D::copySolution( const solution_Type& solution1, solution_Type& solution2 )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::copySolution( solution1, solution2 ) \n";
#endif

    for ( solutionConstIterator_Type i = solution2.begin() ; i != solution2.end() ; ++i )
        if ( solution1.find( i->first ) != solution1.end() )
            *solution2[i->first] = *solution1.find(i->first)->second;
}

void
MultiscaleModel1D::updateBCPhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::updateBCPhysicalSolverVariables() \n";
#endif

    // Update BCInterface solver variables
    M_bc->updatePhysicalSolverVariables();
}

void
MultiscaleModel1D::solve( bc_Type& bc, solution_Type& solution, const std::string& solverType )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::solve() \n";
#endif

    // Re-initialize solution
    copySolution( *M_solution_tn, solution );

    // Subiterate to respect CFL
    UInt subiterationNumber(1);
    Real timeStep = M_data->dataTime()->timeStep();

    Real CFL = M_solver->computeCFL( solution, M_data->dataTime()->timeStep() );
    if ( CFL > M_data->CFLmax() )
    {
        subiterationNumber = std::ceil( CFL / M_data->CFLmax() );
        timeStep /= subiterationNumber;
    }

    if ( M_comm->MyPID() == 0 )
        std::cout << solverType << "  Number of subiterations                  " << subiterationNumber
                                << " ( CFL = " << CFL*timeStep/M_data->dataTime()->timeStep() << " )" << std::endl;

    for ( UInt i(1) ; i <= subiterationNumber ; ++i )
    {
        updateBCPhysicalSolverVariables();
        M_physics->setArea_tn( *solution["A"] );
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
    M_bcPreviousTimeSteps[0][OneDimensional::left][OneDimensional::S]  = M_solver->boundaryValue( solution, OneDimensional::S, OneDimensional::left );
    M_bcPreviousTimeSteps[0][OneDimensional::left][OneDimensional::Q]  = M_solver->boundaryValue( solution, OneDimensional::Q, OneDimensional::left );
    M_bcPreviousTimeSteps[0][OneDimensional::right][OneDimensional::A] = M_solver->boundaryValue( solution, OneDimensional::A, OneDimensional::right );
    M_bcPreviousTimeSteps[0][OneDimensional::right][OneDimensional::S] = M_solver->boundaryValue( solution, OneDimensional::S, OneDimensional::right );
    M_bcPreviousTimeSteps[0][OneDimensional::right][OneDimensional::Q] = M_solver->boundaryValue( solution, OneDimensional::Q, OneDimensional::right );
}

void
MultiscaleModel1D::setupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::setupLinearModel( ) \n";
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
    Debug( 8130 ) << "MultiscaleModel1D::solveLinearModel() \n";
#endif

    if ( !solveLinearSystem )
        return;

    imposePerturbation();

    displayModelStatus( "Solve linear" );
    solve( *M_linearBC, *M_linearSolution, "L1D-" );

    resetPerturbation();

    //This flag avoid recomputation of the same system
    solveLinearSystem = false;
}

void
MultiscaleModel1D::imposePerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::imposePerturbation() \n";
#endif

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            // Find the side to perturb and apply the perturbation
            M_bcDeltaSide = flagConverter( ( *i )->flag( ( *i )->modelGlobalToLocalID( M_ID ) ) );
            M_linearBC->bc( M_bcDeltaSide )->setBCFunction( OneDimensional::first, M_bcBaseDelta );

            // Compute the range
            M_bcDeltaType = M_linearBC->bc( M_bcDeltaSide )->type( OneDimensional::first );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA
            // We replace stress BC with area BC for the perturbed problem
            if ( M_bcDeltaType == OneDimensional::S )
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

            case OneDimensional::S:

                M_bcDelta = M_data->jacobianPerturbationStress();

                break;

            default:

                std::cout << "Warning: bcType \"" << M_bcDeltaType << "\"not available!" << std::endl;

                break;
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
    Debug( 8130 ) << "MultiscaleModel1D::resetPerturbation() \n";
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
    // Previous bc size
    UInt bcPreviousSize( M_bcPreviousTimeSteps.size() );

    // Time container for interpolation
    std::vector< Real > timeContainer( bcPreviousSize, 0 );
    for ( UInt i(0) ; i < bcPreviousSize ; ++i )
        timeContainer[i] = M_globalData->dataTime()->time() - i * M_globalData->dataTime()->timeStep();

    // Lagrange interpolation
    Real bcValue(0);
    Real base(1);

    for ( UInt i(0) ; i < M_bcPreviousTimeSteps.size() ; ++i )
    {
        base = 1;
        for ( UInt j(0) ; j < M_bcPreviousTimeSteps.size() ; ++j )
            if ( j != i )
                base *= (t - timeContainer[j]) / (timeContainer[i] - timeContainer[j]);

        if ( i == 0 )
            bcValue += base * ( M_bcPreviousTimeSteps[i][M_bcDeltaSide][M_bcDeltaType] + M_bcDelta );
        else
            bcValue += base * M_bcPreviousTimeSteps[i][M_bcDeltaSide][M_bcDeltaType];
    }

    return bcValue;
}

#else

Real
MultiscaleModel1D::tangentProblem( const bcSide_Type& bcOutputSide, const bcType_Type& outputType )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MultiscaleModel1D::tangentProblem( bcOutputSide, bcOutputType ) \n";
#endif

    displayModelStatus( "Solve linear" );
    Real jacobianCoefficient(0);

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            // Find the perturbed side
            bcSide_Type bcSide = flagConverter( ( *i )->flag( ( *i )->modelGlobalToLocalID( M_ID ) ) );

            // Perturbation has no effect on the other sides (which also means that dQ/dQ and dP/dP are always zero)
            if ( bcSide != bcOutputSide )
                break;

            // Compute the RHS
            solver_Type::vector_Type rhs( M_feSpace->map() );
            rhs = 0.;

            // Compute the eigenvectors
            data_Type::container2D_Type eigenvalues, leftEigenvector1, leftEigenvector2;
            M_solver->boundaryEigenValuesEigenVectors( bcSide, *M_solution_tn, eigenvalues, leftEigenvector1, leftEigenvector2 );

            UInt bcNode = M_solver->boundaryDOF( bcOutputSide );

            switch ( bcOutputSide )
            {
            case OneDimensional::left:
                switch ( outputType )
                {
                case OneDimensional::Q: // dQ_L/dS_L given by -1 * -1 * ( -L21/L22 )

                    rhs[bcNode] = leftEigenvector2[0] / leftEigenvector2[1];
                    jacobianCoefficient = -solveTangentProblem( rhs, bcNode ) * M_physics->dAdP( M_solver->boundaryValue( *M_solution, OneDimensional::P, bcOutputSide ), M_data->dataTime()->timeStep(), bcNode );

                    break;

                case OneDimensional::S: // dS_L/dQ_L given by -1 * -1 * ( -L22/L21 )

                    jacobianCoefficient = -leftEigenvector2[1] / leftEigenvector2[0]
                                          * M_physics->dPdA( M_solver->boundaryValue( *M_solution, OneDimensional::A, bcOutputSide ), M_data->dataTime()->timeStep(), bcNode );

                    break;

                default:

                    std::cout << "Warning: bcType \"" << outputType << "\"not available!" << std::endl;

                    break;
                }

                break;

           case OneDimensional::right:
                switch ( outputType )
                {
                case OneDimensional::Q: // dQ_R/dS_R given by -1 * -L11/L12

                    rhs[bcNode] = leftEigenvector1[0] / leftEigenvector1[1];
                    jacobianCoefficient = solveTangentProblem( rhs, bcNode ) * M_physics->dAdP( M_solver->boundaryValue( *M_solution, OneDimensional::P, bcOutputSide ), M_data->dataTime()->timeStep(), bcNode );

                    break;

                case OneDimensional::S: // dS_R/dQ_R given by -1 * -L12/L11

                    jacobianCoefficient = leftEigenvector1[1] / leftEigenvector1[0]
                                          * M_physics->dPdA( M_solver->boundaryValue( *M_solution, OneDimensional::A, bcOutputSide ), M_data->dataTime()->timeStep(), bcNode );
                    break;

                default:

                    std::cout << "Warning: bcType \"" << outputType << "\"not available!" << std::endl;

                    break;
                }

                break;

            default:

                std::cout << "Warning: bcSide \"" << bcSide << "\" not available!" << std::endl;

                break;
            }

            // Quit the loop
            break;
        }

    return jacobianCoefficient;
}

Real
MultiscaleModel1D::solveTangentProblem( solver_Type::vector_Type& rhs, const UInt& bcNode )
{
#ifdef HAVE_NEUMANN_VISCOELASTIC_BC
    if ( M_data->viscoelasticWall() )
    {
        // Solve elastic tangent problem
        solver_Type::vector_Type flowRateElasticCorrection( M_feSpace->map() );
        M_linearSolver->solveSystem( rhs, flowRateElasticCorrection, M_solver->massMatrix() );

        // Solve viscoelastic tangent problem
        solver_Type::vector_Type flowRateViscoelasticCorrection( M_feSpace->map() );
        flowRateViscoelasticCorrection = 0.;

        solver_Type::vector_Type flowRateDelta( M_feSpace->map() );
        flowRateDelta = 0.;

        for ( UInt iNode(0); iNode < M_physics->data()->numberOfNodes() ; ++iNode )
        {
            flowRateDelta[iNode] = 1.;
            flowRateViscoelasticCorrection += M_solver->viscoelasticFluxCorrection( *(*M_solution)["A"], flowRateDelta, M_data->dataTime()->timeStep(), *M_bc->handler(), false ) * flowRateElasticCorrection[iNode];
            flowRateDelta[iNode] = 0.;
        }

        return flowRateViscoelasticCorrection[bcNode] + flowRateElasticCorrection[bcNode];
    }
    else
#endif
        return rhs[bcNode];
}

#endif

} // Namespace multiscale
} // Namespace LifeV
