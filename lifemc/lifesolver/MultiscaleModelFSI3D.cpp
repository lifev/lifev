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
 *  @brief File containing the Multiscale Model FSI3D
 *
 *  @date 19-04-2010
 *  @author Paolo Crosetto <paolo.crosetto@epfl.ch>
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleModelFSI3D.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModelFSI3D::MultiscaleModelFSI3D() :
        multiscaleModel_Type           (),
        MultiscaleInterfaceFluid       (),
        M_FSIoperator                  (),
        M_data                         ( new data_Type() ),
        M_exporterFluid                (),
        M_exporterSolid                (),
        M_importerFluid                (),
        M_importerSolid                (),
#ifndef FSI_WITH_EXTERNALPRESSURE
        M_stressCouplingFunction       (),
        M_externalPressureScalar       (),
#endif
        M_externalPressureVector       (),
        M_fluidVelocityAndPressure     (),
        M_fluidDisplacement            (),
        M_solidVelocity                (),
        M_solidDisplacement            (),
        M_fluidVelocityAndPressure_tn  (),
        M_fluidDisplacement_tn         (),
        M_solidVelocity_tn             (),
        M_solidDisplacement_tn         (),
        M_nonLinearRichardsonIteration ( 0 ),
        M_fluidBC                      ( new bcInterface_Type() ),
        M_solidBC                      ( new bcInterface_Type() ),
        M_harmonicExtensionBC          ( new bcInterface_Type() ),
        M_linearBC                     (),
        M_linearRHS                    (),
        M_linearSolution               ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::MultiscaleModelFSI3D() \n";
#endif

    M_type = FSI3D;

    BlockPrecFactory::instance().registerProduct("AdditiveSchwarz",   &MonolithicBlockMatrix::createAdditiveSchwarz) ;
    BlockPrecFactory::instance().registerProduct("AdditiveSchwarzRN", &MonolithicBlockMatrixRN::createAdditiveSchwarzRN ) ;
    BlockPrecFactory::instance().registerProduct("ComposedDN",        &MonolithicBlockComposedDN::createComposedDN ) ;
    BlockPrecFactory::instance().registerProduct("ComposedDN2",       &MonolithicBlockComposedDN::createComposedDN2 );
    BlockPrecFactory::instance().registerProduct("ComposedDNND",      &MonolithicBlockComposedDNND::createComposedDNND);

    MonolithicBlockMatrix::Factory_Type::instance().registerProduct("AdditiveSchwarz",   &MonolithicBlockMatrix::createAdditiveSchwarz ) ;
    MonolithicBlockMatrix::Factory_Type::instance().registerProduct("AdditiveSchwarzRN", &MonolithicBlockMatrixRN::createAdditiveSchwarzRN ) ;

    FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "linearVenantKirchhof", &FSIOperator::createVenantKirchhoffLinear );
    FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "exponential", &FSIOperator::createExponentialMaterialNonLinear );
    FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "neoHookean", &FSIOperator::createNeoHookeanMaterialNonLinear );
    FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct( "nonLinearVenantKirchhof", &FSIOperator::createVenantKirchhoffNonLinear );


    //FSIOperator_Type::solid_Type::StructureSolverFactory::instance().registerProduct( "linearVenantKirchhof", &FSIOperator_Type::createLinearStructure );
//    FSIOperator_Type::solid_Type::StructureSolverFactory_Type::instance().registerProduct( "nonLinearVenantKirchhof", &FSIOperator_Type::createNonLinearStructure );
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModelFSI3D::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::setupData( fileName ) \n";
#endif

    multiscaleModel_Type::setupData( fileName );

    GetPot dataFile( fileName );

    // Load data
    M_data->setup( dataFile );
    if ( M_globalData.get() )
        setupGlobalData( fileName );

    // Create FSI
    M_FSIoperator = FSIOperatorPtr_Type( FSIOperator::FSIFactory_Type::instance().createObject( M_data->method() ) );

    // Setup Communicator
    setupCommunicator();

    // Set data
    M_FSIoperator->setData( M_data );
    M_FSIoperator->setDataFile( dataFile );

    // Setup Boundary Conditions from file
    setupBC( fileName );

    // Setup exporter
    if ( M_FSIoperator->isFluid() )
    {
        setupExporter( M_exporterFluid, dataFile );
        setupImporter( M_importerFluid, dataFile );
    }
    if ( M_FSIoperator->isSolid() )
    {
        setupExporter( M_exporterSolid, dataFile, "_Solid" );
        setupImporter( M_importerSolid, dataFile, "_Solid" );
    }
}

void
MultiscaleModelFSI3D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::setupProblem() \n";
#endif

    // Mesh transformation (before partitioning, ideally should be done after for scalability)
    //M_FSIoperator->fluidMesh().transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );
    M_FSIoperator->solidMesh().meshTransformer().transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    // Mesh partitioning
    M_FSIoperator->partitionMeshes();

    // Mesh transformation (after partitioning - not working for solid)
    M_FSIoperator->fluidMeshPart().meshPartition()->meshTransformer().transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );
    //M_FSIoperator->solidMeshPart().meshPartition()->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    // Setup FEspace & DOF
    M_FSIoperator->setupFEspace();
    M_FSIoperator->setupDOF();

    // Setup FSI Interface Boundary Conditions (by giving the operator to BCInterface)
    M_fluidBC->setPhysicalSolver( M_FSIoperator );
    M_solidBC->setPhysicalSolver( M_FSIoperator );
    M_harmonicExtensionBC->setPhysicalSolver( M_FSIoperator );

    // Setup Fluid & Solid solver
    M_FSIoperator->setupFluidSolid( M_FSIoperator->imposedFluxes() );
    M_FSIoperator->setupSystem();

    // Setup Exporters
    if ( M_FSIoperator->isFluid() )
        setExporterFluid( M_exporterFluid );

    if ( M_FSIoperator->isSolid() )
        setExporterSolid( M_exporterSolid );

    //Setup solution
    initializeSolution();

    //Setup linear model
    setupLinearModel();

    M_FSIoperator->fluid().setupPostProc(); //this has to be called if we want to initialize the postProcess
}

void
MultiscaleModelFSI3D::buildModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::buildModel() \n";
#endif

    // Display data
//    if ( M_comm->MyPID() == 0 )
//        M_data->showMe();

    // Update BCInterface solver variables
    updateBC();

    M_FSIoperator->buildSystem();
    M_FSIoperator->updateSystem();
}

void
MultiscaleModelFSI3D::updateModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::updateModel() \n";
#endif

    // Update BCInterface solver variables
    updateBC();

    // Update System
    M_FSIoperator->updateSystem();

    // M_FSIoperator->fluidTimeAdvance()->shiftRight(*M_fluidVelocityAndPressure);
    // M_FSIoperator->solidTimeAdvance()->shiftRight(*M_solidDisplacement);
    //this is already done by updateSolution!
    M_FSIoperator->updateSolution( *M_stateVariable );
    //calls the shiftRight

    // Update solution at time n
    // *M_fluidVelocityAndPressure_tn = *M_fluidVelocityAndPressure;
    // *M_fluidDisplacement_tn        = *M_fluidDisplacement;
    // *M_solidVelocity_tn            = *M_solidVelocity;
    // *M_solidDisplacement_tn        = *M_solidDisplacement;


    // Parameters for Multiscale subiterations
    boost::dynamic_pointer_cast< FSIMonolithic > ( M_FSIoperator )->precPtrView()->setRecompute( 1, true );
    M_nonLinearRichardsonIteration = 0;
}

void
MultiscaleModelFSI3D::solveModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::solveModel() \n";
#endif

    displayModelStatus( "Solve" );

    // Initialize all the quantities in the solver to time tn
// <<<<<<< HEAD
    //    M_FSIoperator->initialize( M_fluidVelocityAndPressure_tn, M_fluidDisplacement_tn, M_solidVelocity_tn, M_solidDisplacement_tn );
// #else
     // vectorPtr_Type fluidVelocityAndPressureWithoutExternalPressure_tn( new vector_Type( *M_fluidVelocityAndPressure_tn - *M_externalPressureVector ) );
     // std::vector<vectorPtr_Type> vel(1), fluidDisp(1), solidVel(1), solidDisp(1);
     // vel[0]=fluidVelocityAndPressureWithoutExternalPressure_tn;
     // fluidDisp[0]=M_fluidDisplacement_tn;
     // solidVel[0]=M_solidVelocity_tn;
     // solidDisp[0]=M_solidDisplacement_tn;
     // M_FSIoperator->initialize(veli, fluidDisp, solidDisp);
// #endif
// =======
// #ifdef FSI_WITH_EXTERNALPRESSURE
//     M_FSIoperator->initialize( M_fluidVelocityAndPressure_tn, M_fluidDisplacement_tn, M_solidVelocity_tn, M_solidDisplacement_tn );
// #else
//     vectorPtr_Type fluidVelocityAndPressureWithoutExternalPressure_tn( new vector_Type( *M_fluidVelocityAndPressure_tn - *M_externalPressureVector ) );
//     M_FSIoperator->initialize( fluidVelocityAndPressureWithoutExternalPressure_tn, M_fluidDisplacement_tn, M_solidVelocity_tn, M_solidDisplacement_tn );
// #endif
// >>>>>>> master

    if ( M_nonLinearRichardsonIteration != 0 )
    {
        M_FSIoperator->resetRHS();
        M_FSIoperator->updateRHS();
        M_FSIoperator->applyBoundaryConditions();
    }

    // Non-linear Richardson solver
    UInt maxSubIterationNumber = M_data->maxSubIterationNumber();

// #ifdef FSI_WITH_EXTERNALPRESSURE
//     vectorPtr_Type solution( new vector_Type( *M_fluidVelocityAndPressure_tn ) );
// #else
//     vectorPtr_Type solution( new vector_Type( *fluidVelocityAndPressureWithoutExternalPressure_tn ) );
// #endif
    vectorPtr_Type solution( new vector_Type( *M_FSIoperator->fluidTimeAdvance()->stencil()[0] ) );// solution at the previous time step
    //why not an extrapolation?

    solution->spy("sol");

    std::ofstream outRes; // Unuseful variable - NonLinearRichardson interface must be clean..

    NonLinearRichardson( *solution, *M_FSIoperator,
                          M_data->absoluteTolerance(), M_data->relativeTolerance(),
                          maxSubIterationNumber, M_data->errorTolerance(),
                          M_data->NonLinearLineSearch(),
                          outRes, M_data->dataFluid()->dataTime()->time(),
                          M_nonLinearRichardsonIteration
                       );

    //I guess this should go to the updateModel method
    //M_FSIoperator->updateSolution( *solution );
    M_stateVariable.reset( new vector_Type(*solution));

    // Parameters for Multiscale subiterations
    boost::dynamic_pointer_cast< FSIMonolithic > ( M_FSIoperator )->precPtrView()->setRecompute( 1, false );
    M_nonLinearRichardsonIteration = 1;
}

void
MultiscaleModelFSI3D::updateSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::updateSolution() \n";
#endif

// TODO Fix postprocessing (remove the 4 variables and use M_stateVariables)
// TODO Remove M_solidVelocity
    if ( M_FSIoperator->isFluid() )
    {
        *M_fluidDisplacement      = M_FSIoperator->meshDisp();
        //M_FSIoperator->exportFluidDisplacement( *M_fluidDisplacement );
        M_FSIoperator->exportFluidVelocityAndPressure( *M_fluidVelocityAndPressure );
        #ifndef FSI_WITH_EXTERNALPRESSURE
        *M_fluidVelocityAndPressure += *M_externalPressureVector;
        #endif
    }

    if ( M_FSIoperator->isSolid() )
    {
        M_FSIoperator->exportSolidDisplacement( *M_solidDisplacement );
        M_FSIoperator->exportSolidVelocity( *M_solidVelocity );
    }
}

void
MultiscaleModelFSI3D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::saveSolution() \n";
#endif

    if ( M_FSIoperator->isFluid() )
        M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->time() );

    if ( M_FSIoperator->isSolid() )
        M_exporterSolid->postProcess( M_data->dataSolid()->dataTime()->time() );

#ifdef HAVE_HDF5
    if ( M_data->dataFluid()->dataTime()->isLastTimeStep() )
    {
        if ( M_FSIoperator->isFluid() )
            M_exporterFluid->closeFile();
        if ( M_FSIoperator->isSolid() )
            M_exporterSolid->closeFile();
    }
#endif

}

void
MultiscaleModelFSI3D::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        multiscaleModel_Type::showMe();

        std::cout << "FSI method          = " << M_data->method() << std::endl << std::endl;

        std::cout << "Velocity FE order   = " << M_FSIoperator->dataFluid()->uOrder() << std::endl
                  << "Pressure FE order   = " << M_FSIoperator->dataFluid()->pOrder() << std::endl
                  << "Structure FE order  = " << M_FSIoperator->dataSolid()->order() << std::endl<< std::endl;

        std::cout << "Velocity DOF        = " << 3 * M_FSIoperator->uFESpace().dof().numTotalDof() << std::endl
                  << "Pressure DOF        = " << M_FSIoperator->pFESpace().dof().numTotalDof() << std::endl
                  << "Harmonic ext. DOF   = " << M_FSIoperator->mmFESpace().dof().numTotalDof() << std::endl
                  << "Structure DOF       = " << M_FSIoperator->dFESpace().dof().numTotalDof() << std::endl << std::endl;

        std::cout << "Fluid mesh maxH     = " << MeshUtility::MeshStatistics::computeSize( *( M_FSIoperator->uFESpace().mesh() ) ).maxH << std::endl
                  << "Fluid mesh meanH    = " << MeshUtility::MeshStatistics::computeSize( *( M_FSIoperator->uFESpace().mesh() ) ).meanH << std::endl
                  << "Solid mesh maxH     = " << MeshUtility::MeshStatistics::computeSize( *( M_FSIoperator->dFESpace().mesh() ) ).maxH << std::endl
                  << "Solid mesh meanH    = " << MeshUtility::MeshStatistics::computeSize( *( M_FSIoperator->dFESpace().mesh() ) ).meanH << std::endl << std::endl;
    }
}

Real
MultiscaleModelFSI3D::checkSolution() const
{
    return M_fluidDisplacement->norm2() + M_fluidVelocityAndPressure->norm2() + M_solidDisplacement->norm2() + M_solidVelocity->norm2();
}



// ===================================================
// MultiscaleInterfaceFluid Methods
// ===================================================
void
MultiscaleModelFSI3D::imposeBoundaryFlowRate( const bcFlag_Type& flag, const function_Type& function )
{
    BCFunctionBase base;
    base.setFunction( function );

    M_fluidBC->handler()->addBC( "CouplingFlowRate_Model_" + number2string( M_ID ) + "_Flag_" + number2string( flag ), flag, Flux, Full, base, 3 );
}

void
MultiscaleModelFSI3D::imposeBoundaryStress( const bcFlag_Type& flag, const function_Type& function )
{
    BCFunctionBase base;
#ifdef FSI_WITH_EXTERNALPRESSURE
    base.setFunction( function );
#else
    FSI3DCouplingFunction couplingFunction( function, M_data->dataSolid()->externalPressure() );
    M_stressCouplingFunction.push_back( couplingFunction );
    base.setFunction( boost::bind( &FSI3DCouplingFunction::function, M_stressCouplingFunction.back(), _1, _2, _3, _4, _5 ) );
#endif
    M_fluidBC->handler()->addBC( "CouplingStress_Model_" + number2string( M_ID ) + "_Flag_" + number2string( flag ), flag, Natural, Normal, base );
}

Real
MultiscaleModelFSI3D::boundaryDeltaFlowRate( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    solveLinearModel( solveLinearSystem );

    return M_FSIoperator->fluid().flux( flag, *M_linearSolution );
}

Real
MultiscaleModelFSI3D::boundaryDeltaStress( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    solveLinearModel( solveLinearSystem );

    if ( M_linearBC->findBCWithFlag( flag ).type() == Flux )
        return -M_FSIoperator->fluid().lagrangeMultiplier( flag, *M_linearBC, *M_linearSolution );
    else
        return -M_FSIoperator->fluid().pressure( flag, *M_linearSolution );
}

// ===================================================
// Get Methods
// ===================================================
Real
MultiscaleModelFSI3D::boundaryPressure( const bcFlag_Type& flag ) const
{
    if ( M_fluidBC->handler()->findBCWithFlag( flag ).type() == Flux )
#ifdef FSI_WITH_EXTERNALPRESSURE
        return M_FSIoperator->fluid().lagrangeMultiplier( flag, *M_fluidBC->handler(), M_FSIoperator->solution() );
#else
        return M_FSIoperator->fluid().lagrangeMultiplier( flag, *M_fluidBC->handler(), M_FSIoperator->solution() ) + M_externalPressureScalar;
#endif
    else
#ifdef FSI_WITH_EXTERNALPRESSURE
        return M_FSIoperator->fluid().pressure( flag, M_FSIoperator->solution() );
#else
        return M_FSIoperator->fluid().pressure( flag, M_FSIoperator->solution() ) + M_externalPressureScalar;
#endif
}

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleModelFSI3D::setupGlobalData( const std::string& fileName )
{
    GetPot dataFile( fileName );

    // Global data time
    M_data->dataFluid()->setTimeData( M_globalData->dataTime() );
    M_data->dataSolid()->setTimeData( M_globalData->dataTime() );

    // Fluid global physical quantities
    if ( !dataFile.checkVariable( "fluid/physics/density" ) )
        M_data->dataFluid()->setDensity( M_globalData->fluidDensity() );
    if ( !dataFile.checkVariable( "fluid/physics/viscosity" ) )
        M_data->dataFluid()->setViscosity( M_globalData->fluidViscosity() );

    // Solid global physical quantities
    if ( !dataFile.checkVariable( "solid/physics/density" ) )
        M_data->dataSolid()->setDensity( M_globalData->solidDensity() );
    if ( !dataFile.checkVariable( "solid/physics/externalPressure" ) )
        M_data->dataSolid()->setExternalPressure( M_globalData->solidExternalPressure() );

    std::vector< UInt > materialFlags;
    if ( !dataFile.checkVariable( "solid/physics/material_flag" ) )
        materialFlags.push_back( 1 );
    else
        for ( UInt i( 0 ) ; i < dataFile.vector_variable_size( "solid/physics/material_flag" ) ; ++i )
            materialFlags.push_back( dataFile( "solid/physics/material_flag", 1, i ) );

    for ( std::vector< UInt >::const_iterator i = materialFlags.begin(); i != materialFlags.end() ; ++i )
    {
        if ( !dataFile.checkVariable( "solid/physics/poisson" ) )
            M_data->dataSolid()->setPoisson( M_globalData->solidPoissonCoefficient(), *i );
        if ( !dataFile.checkVariable( "solid/physics/young" ) )
            M_data->dataSolid()->setYoung( M_globalData->solidYoungModulus(), *i );
    }
}

void
MultiscaleModelFSI3D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::initializeSolution() \n";
#endif

    // Initialize the external pressure vector
    vector_Type fluidPressure( M_FSIoperator->pFESpace().mapPtr(), Unique );
    fluidPressure = M_data->dataSolid()->externalPressure();

    M_externalPressureVector.reset( new vector_Type( *M_fluidVelocityAndPressure, Unique, Zero ) );
    M_externalPressureVector->subset( fluidPressure, fluidPressure.map(), static_cast <UInt> ( 0 ), static_cast <UInt>( 3 * M_FSIoperator->uFESpace().dof().numTotalDof() ) );

#ifndef FSI_WITH_EXTERNALPRESSURE
    // Initialize the external pressure scalar
    M_externalPressureScalar = M_data->dataSolid()->externalPressure();
    M_data->dataSolid()->setExternalPressure( 0 );
    updateBC();
#endif

    // Initialize solution at time tn
    M_fluidVelocityAndPressure_tn.resize(  M_FSIoperator->fluidTimeAdvance()->size()  );
    M_fluidDisplacement_tn.resize( M_FSIoperator->ALETimeAdvance()->size() ) ;
    //M_solidVelocity_tn.reset( new vector_Type( *M_solidVelocity ) );
    M_solidDisplacement_tn.resize( M_FSIoperator->solidTimeAdvance()->size() );
    vectorPtr_Type firstFluidDisp( new vector_Type(*M_FSIoperator->couplingVariableMap()) );

    vectorPtr_Type solution( new vector_Type(*M_FSIoperator->couplingVariableMap()) );

    if ( multiscaleProblemStep > 0 )
    {
        M_importerFluid->setMeshProcId( M_FSIoperator->uFESpace().mesh(), M_FSIoperator->uFESpace().map().comm().MyPID() );
        M_importerSolid->setMeshProcId( M_FSIoperator->dFESpace().mesh(), M_FSIoperator->dFESpace().map().comm().MyPID() );

        M_importerFluid->addVariable( IOData_Type::VectorField, "Displacement (fluid)", M_FSIoperator->mmFESpacePtr(), M_fluidDisplacement, static_cast <UInt> (0) );
        M_importerFluid->addVariable( IOData_Type::VectorField, "Velocity (fluid)",     M_FSIoperator->uFESpacePtr(),  M_fluidVelocityAndPressure, static_cast <UInt> (0) );
        M_importerFluid->addVariable( IOData_Type::ScalarField, "Pressure (fluid)",     M_FSIoperator->pFESpacePtr(),  M_fluidVelocityAndPressure, static_cast <UInt> (3 * M_FSIoperator->uFESpace().dof().numTotalDof() ) );

        M_importerSolid->addVariable( IOData_Type::VectorField, "Displacement (solid)", M_FSIoperator->dFESpacePtr(),  M_solidDisplacement, static_cast <UInt> (0) );
        M_importerSolid->addVariable( IOData_Type::VectorField, "Velocity (solid)",     M_FSIoperator->dFESpacePtr(),  M_solidVelocity,     static_cast <UInt> (0) );

        vectorPtr_Type temporaryVector;


        for(UInt iteration(0); iteration<M_FSIoperator->fluidTimeAdvance()->size() ; ++iteration)
        {
            // Import
            M_exporterFluid->setTimeIndex( M_importerFluid->importFromTime( M_data->dataFluid()->dataTime()->initialTime()-iteration*M_data->dataFluid()->dataTime()->timeStep() ) + 1 );

#ifdef HAVE_HDF5
        if ( M_FSIoperator->isFluid() )
            M_importerFluid->closeFile();
        if ( M_FSIoperator->isSolid() )
            M_importerSolid->closeFile();
#endif

            // Assemble the Monolithic solution
            // Add fluid
            temporaryVector.reset( new vector_Type( *M_fluidVelocityAndPressure, Unique, Zero ));
            *solution = *temporaryVector;

            M_fluidVelocityAndPressure_tn[iteration]=solution;

        }



        for(UInt iteration(0); iteration<M_FSIoperator->solidTimeAdvance()->size() ; ++iteration)
        {
            M_exporterSolid->setTimeIndex( M_importerSolid->importFromTime( M_data->dataSolid()->dataTime()->initialTime()-iteration*M_data->dataSolid()->dataTime()->timeStep() ) + 1 );

            // Add solid
            UInt offset = boost::dynamic_pointer_cast< FSIMonolithic > ( M_FSIoperator )->offset();

            // Assemble the Monolithic solution
            temporaryVector.reset( new vector_Type(*M_FSIoperator->couplingVariableMap(), Unique, Zero));
            temporaryVector->subset( *M_solidDisplacement, M_solidDisplacement->map(), static_cast<UInt> ( 0 ), offset );
            *temporaryVector /= M_FSIoperator->solid().rescaleFactor();

            *solution = *temporaryVector;
            M_solidDisplacement_tn[iteration]=solution;
        }

        Int convectiveTerm=0; //loadInitSol;
        if(!M_data->dataFluid()->domainVelImplicit())
            convectiveTerm=1;

        for(UInt iteration(0); iteration<M_FSIoperator->ALETimeAdvance()->size()+convectiveTerm+1 ; ++iteration)
        {
            temporaryVector.reset( new vector_Type(*M_FSIoperator->couplingVariableMap(), Unique, Zero));
            M_exporterFluid->setTimeIndex( M_importerFluid->importFromTime( M_data->dataFluid()->dataTime()->initialTime()-iteration*M_data->dataFluid()->dataTime()->timeStep() ) + 1 );
            // Add harmonic extension
            if ( !M_data->method().compare("monolithicGI") )
            {
                UInt offset = boost::dynamic_pointer_cast< FSIMonolithicGI > ( M_FSIoperator )->mapWithoutMesh().map( Unique )->NumGlobalElements();

                if(iteration == 0)
                {
                    temporaryVector->subset(*M_fluidDisplacement, M_fluidDisplacement->map(), (UInt)0, offset);
                    if(!convectiveTerm)
                    {
                        *firstFluidDisp = *M_fluidDisplacement;
                    }
                }
                else
                    if(convectiveTerm && iteration == 1)
                        *firstFluidDisp = *M_fluidDisplacement;
                    else
                        M_fluidDisplacement_tn[iteration]=(M_fluidDisplacement);
            }
        }
        *M_fluidVelocityAndPressure_tn[0]+=*M_solidDisplacement_tn[0];
        *M_fluidVelocityAndPressure_tn[0]+=*temporaryVector;
    }
    else
    {
        for(UInt iteration(0); iteration<M_FSIoperator->fluidTimeAdvance()->size() ; ++iteration)
        {
            // Initialize solution
            *solution = *M_externalPressureVector;
            // Initialize solution at time tn
            M_fluidVelocityAndPressure_tn[iteration]=solution;
        }
        for(UInt iteration(0); iteration<M_FSIoperator->ALETimeAdvance()->size() ; ++iteration)
        {
            *M_fluidDisplacement        = 0.0;
            M_fluidDisplacement_tn[iteration]=M_fluidDisplacement;
        }
        for(UInt iteration(0); iteration<M_FSIoperator->solidTimeAdvance()->size() ; ++iteration)
        {
            *solution                        = 0.0;
            M_solidDisplacement_tn[iteration]=solution;
        }
    }

    vectorPtr_Type fluidVelocityAndPressureWithoutExternalPressure_tn( new vector_Type( *M_fluidVelocityAndPressure_tn[0] - *M_externalPressureVector ) );
    //    std::vector<vectorPtr_Type> vel(1), fluidDisp(1), solidDisp(1);
    //    fluidDisp[0]=M_fluidDisplacement_tn;
    //    solidDisp[0]=M_solidDisplacement_tn;
    M_FSIoperator->initializeTimeAdvance(M_fluidVelocityAndPressure_tn, M_solidDisplacement_tn, M_fluidDisplacement_tn);

    if(!M_data->dataFluid()->domainVelImplicit() &&  multiscaleProblemStep == 0)
      {
          //The following is needed because (and if) of the extrapolation of the fluid domain velocity is used, i.e. M_domainVelImplicit==false
          M_FSIoperator->ALETimeAdvance()->updateRHSFirstDerivative( M_data->dataSolid()->dataTime()->timeStep() );
          M_FSIoperator->ALETimeAdvance()->shiftRight(*firstFluidDisp);
      }

    // Initialize all the quantities in the solver to time tn

    //to Cristiano: Check below
// #ifdef FSI_WITH_EXTERNALPRESSURE
//     *M_fluidVelocityAndPressure_tn[0]=fluidVelocityAndPressure_tn;
//     M_FSIoperator->initialize( vel, fluidDisp, solidDisp );
//     //M_FSIoperator->initialize( M_fluidVelocityAndPressure_tn, M_fluidDisplacement_tn, M_solidVelocity_tn, M_solidDisplacement_tn );
// #else
//     vel[0]=fluidVelocityAndPressureWithoutExternalPressure_tn;
//     M_FSIoperator->initialize( vel, fluidDisp, solidDisp );
// #endif
}

void
MultiscaleModelFSI3D::setupCommunicator()
{
    M_FSIoperator->setFluid( true );
    M_FSIoperator->setSolid( true );

    M_FSIoperator->setFluidLeader( 0 );
    M_FSIoperator->setSolidLeader( 0 );

    M_FSIoperator->setComm( M_comm, M_comm );
}

void
MultiscaleModelFSI3D::setupBC( const std::string& fileName )
{
    if ( M_FSIoperator->isFluid() )
    {
        M_fluidBC->createHandler();
        M_fluidBC->fillHandler( fileName, "fluid" );

        M_FSIoperator->setFluidBC( M_fluidBC->handler() );

        M_harmonicExtensionBC->createHandler();
        M_harmonicExtensionBC->fillHandler( fileName, "mesh_motion" );

        M_FSIoperator->setHarmonicExtensionBC( M_harmonicExtensionBC->handler() );
    }

    if ( M_FSIoperator->isSolid() )
    {
        M_solidBC->createHandler();
        M_solidBC->fillHandler( fileName, "solid" );

        M_FSIoperator->setSolidBC( M_solidBC->handler() );
    }
}

void
MultiscaleModelFSI3D::updateBC()
{
    if ( M_FSIoperator->isFluid() )
    {
        M_fluidBC->updatePhysicalSolverVariables();
        M_harmonicExtensionBC->updatePhysicalSolverVariables();
    }
    if ( M_FSIoperator->isSolid() )
    {
        M_solidBC->updatePhysicalSolverVariables();
    }
}

void
MultiscaleModelFSI3D::setupExporter( IOFilePtr_Type& exporter, const GetPot& dataFile, const std::string& label )
{
    const std::string exporterType = dataFile( "exporter/type", "ensight" );
#ifdef HAVE_HDF5
    if ( exporterType.compare( "hdf5" ) == 0 )
        exporter.reset( new hdf5IOFile_Type() );
    else
#endif
        exporter.reset( new ensightIOFile_Type() );

    exporter->setDataFromGetPot( dataFile );
    exporter->setPrefix( multiscaleProblemPrefix + "_Model_" + number2string( M_ID ) +  label + "_" + number2string( multiscaleProblemStep ) );
    exporter->setPostDir( multiscaleProblemFolder );
}

void
MultiscaleModelFSI3D::setupImporter( IOFilePtr_Type& importer, const GetPot& dataFile, const std::string& label )
{
    const std::string importerType = dataFile( "exporter/type", "ensight" );
#ifdef HAVE_HDF5
    if ( importerType.compare( "hdf5" ) == 0 )
        importer.reset( new hdf5IOFile_Type() );
    else
#endif
        importer.reset( new ensightIOFile_Type() );

    importer->setDataFromGetPot( dataFile );
    importer->setPrefix( multiscaleProblemPrefix + "_Model_" + number2string( M_ID ) +  label + "_" + number2string( multiscaleProblemStep - 1 ) );
    importer->setPostDir( multiscaleProblemFolder );
}

void
MultiscaleModelFSI3D::setExporterFluid( const IOFilePtr_Type& exporter )
{
    M_fluidVelocityAndPressure.reset( new vector_Type( M_FSIoperator->fluid().getMap(),  M_exporterFluid->mapType() ) );
    M_fluidDisplacement.reset       ( new vector_Type( M_FSIoperator->mmFESpace().map(), M_exporterFluid->mapType() ) );

    exporter->setMeshProcId( M_FSIoperator->uFESpace().mesh(), M_FSIoperator->uFESpace().map().comm().MyPID() );

    exporter->addVariable( IOData_Type::VectorField, "Displacement (fluid)", M_FSIoperator->mmFESpacePtr(), M_fluidDisplacement,        static_cast<UInt> ( 0 ) );
    exporter->addVariable( IOData_Type::VectorField, "Velocity (fluid)",     M_FSIoperator->uFESpacePtr(),  M_fluidVelocityAndPressure, static_cast<UInt> ( 0 ) );
    exporter->addVariable( IOData_Type::ScalarField, "Pressure (fluid)",     M_FSIoperator->pFESpacePtr(),  M_fluidVelocityAndPressure, static_cast<UInt> (3 * M_FSIoperator->uFESpace().dof().numTotalDof() ) );
}

void
MultiscaleModelFSI3D::setExporterSolid( const IOFilePtr_Type& exporter )
{
    M_solidDisplacement.reset( new vector_Type( M_FSIoperator->dFESpace().map(), M_exporterSolid->mapType() ) );
    M_solidVelocity.reset    ( new vector_Type( M_FSIoperator->dFESpace().map(), M_exporterSolid->mapType() ) );

    exporter->setMeshProcId( M_FSIoperator->dFESpace().mesh(), M_FSIoperator->dFESpace().map().comm().MyPID() );

    exporter->addVariable( IOData_Type::VectorField, "Displacement (solid)", M_FSIoperator->dFESpacePtr(), M_solidDisplacement, static_cast<UInt> ( 0 ) );
    exporter->addVariable( IOData_Type::VectorField, "Velocity (solid)",     M_FSIoperator->dFESpacePtr(), M_solidVelocity,     static_cast<UInt> ( 0 ) );
}

void
MultiscaleModelFSI3D::setupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::setupLinearModel() \n";
#endif

    // The linear BCHandler is a copy of the original BCHandler with all BCFunctions giving zero
    M_linearBC.reset( new bc_Type( *M_fluidBC->handler() ) );

    // Set all the BCFunctions to zero
    BCFunctionBase bcBaseDeltaZero;
    bcBaseDeltaZero.setFunction( boost::bind( &MultiscaleModelFSI3D::bcFunctionDeltaZero, this, _1, _2, _3, _4, _5 ) );

    for ( bc_Type::bcBaseIterator_Type i = M_linearBC->begin() ; i != M_linearBC->end() ; ++i )
        i->setBCFunction( bcBaseDeltaZero );

    // Setup linear solution & the RHS
    M_linearSolution.reset( new vector_Type( *M_FSIoperator->couplingVariableMap() ) );
    M_linearRHS.reset( new vector_Type( *M_FSIoperator->couplingVariableMap() ) );
}

void
MultiscaleModelFSI3D::updateLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::updateLinearModel() \n";
#endif

    //Create the RHS
    *M_linearRHS *= 0;
    M_FSIoperator->bcManageVectorRHS( M_linearBC, *M_linearRHS );
}

void
MultiscaleModelFSI3D::solveLinearModel( bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::solveLinearModel() \n";
#endif

    if ( !solveLinearSystem )
        return;

    imposePerturbation();

    updateLinearModel();

    //Solve the linear problem
    displayModelStatus( "Solve linear" );
    M_FSIoperator->solveJac( *M_linearSolution, *M_linearRHS, 0. );

    resetPerturbation();

    //This flag avoid recomputation of the same system
    solveLinearSystem = false;
}

void
MultiscaleModelFSI3D::imposePerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::imposePerturbation() \n";
#endif

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            BCFunctionBase bcBaseDeltaOne;
            bcBaseDeltaOne.setFunction( boost::bind( &MultiscaleModelFSI3D::bcFunctionDeltaOne, this, _1, _2, _3, _4, _5 ) );

            M_linearBC->findBCWithFlag( ( *i )->flag( ( *i )->modelGlobalToLocalID( M_ID ) ) ).setBCFunction( bcBaseDeltaOne );

            break;
        }
}

void
MultiscaleModelFSI3D::resetPerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::resetPerturbation() \n";
#endif

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            BCFunctionBase bcBaseDeltaZero;
            bcBaseDeltaZero.setFunction( boost::bind( &MultiscaleModelFSI3D::bcFunctionDeltaZero, this, _1, _2, _3, _4, _5 ) );

            M_linearBC->findBCWithFlag( ( *i )->flag( ( *i )->modelGlobalToLocalID( M_ID ) ) ).setBCFunction( bcBaseDeltaZero );

            break;
        }
}

} // Namespace multiscale
} // Namespace LifeV
