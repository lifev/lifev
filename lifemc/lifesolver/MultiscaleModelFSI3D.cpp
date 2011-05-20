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

    FSIOperator_Type::solid_Type::StructureSolverFactory::instance().registerProduct( "linearVenantKirchhof", &FSIOperator_Type::createLinearStructure );
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
    M_FSIoperator->solidMesh().transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    // Mesh partitioning
    M_FSIoperator->partitionMeshes();

    // Mesh transformation (after partitioning - not working for solid)
    M_FSIoperator->fluidMeshPart().meshPartition()->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );
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

    // Update solution at time n
    *M_fluidVelocityAndPressure_tn = *M_fluidVelocityAndPressure;
    *M_fluidDisplacement_tn        = *M_fluidDisplacement;
    *M_solidVelocity_tn            = *M_solidVelocity;
    *M_solidDisplacement_tn        = *M_solidDisplacement;

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
    M_FSIoperator->initialize( M_fluidVelocityAndPressure_tn, M_fluidDisplacement_tn, M_solidVelocity_tn, M_solidDisplacement_tn );

    if ( M_nonLinearRichardsonIteration != 0 )
    {
        M_FSIoperator->resetRHS();
        M_FSIoperator->updateRHS();
        M_FSIoperator->applyBoundaryConditions();
    }

    // Non-linear Richardson solver
    UInt maxSubIterationNumber = M_data->maxSubIterationNumber();
    vectorPtr_Type solution( new vector_Type( *M_fluidVelocityAndPressure_tn ) );
    std::ofstream outRes; // Unuseful variable - NonLinearRichardson interface must be clean..

    NonLinearRichardson( *solution, *M_FSIoperator,
                          M_data->absoluteTolerance(), M_data->relativeTolerance(),
                          maxSubIterationNumber, M_data->errorTolerance(),
                          M_data->NonLinearLineSearch(),
                          outRes, M_data->dataFluid()->dataTime()->time(),
                          M_nonLinearRichardsonIteration
                       );

    M_FSIoperator->updateSolution( *solution );

    // Parameters for Multiscale subiterations
    boost::dynamic_pointer_cast< FSIMonolithic > ( M_FSIoperator )->precPtrView()->setRecompute( 1, false );
    M_nonLinearRichardsonIteration = 1;
}

void
MultiscaleModelFSI3D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::saveSolution() \n";
#endif

    updateSolution();

    if ( M_FSIoperator->isFluid() )
        M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->time() );

    if ( M_FSIoperator->isSolid() )
        M_exporterSolid->postProcess( M_data->dataSolid()->getdataTime()->time() );

#ifdef HAVE_HDF5
    if ( M_data->dataFluid()->dataTime()->isLastTimeStep() )
    {
        if ( M_FSIoperator->isFluid() )
            ( multiscaleDynamicCast< hdf5IOFile_Type >( M_exporterFluid ) )->closeFile();
        if ( M_FSIoperator->isSolid() )
            ( multiscaleDynamicCast< hdf5IOFile_Type >( M_exporterSolid ) )->closeFile();
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
                  << "Structure FE order  = " << M_FSIoperator->dataSolid()->getOrder() << std::endl<< std::endl;

        std::cout << "Velocity DOF        = " << 3 * M_FSIoperator->uFESpace().dof().numTotalDof() << std::endl
                  << "Pressure DOF        = " << M_FSIoperator->pFESpace().dof().numTotalDof() << std::endl
                  << "Harmonic ext. DOF   = " << M_FSIoperator->mmFESpace().dof().numTotalDof() << std::endl
                  << "Structure DOF       = " << M_FSIoperator->dFESpace().dof().numTotalDof() << std::endl << std::endl;

        std::cout << "Fluid mesh maxH     = " << M_FSIoperator->uFESpace().mesh()->maxH() << std::endl
                  << "Fluid mesh meanH    = " << M_FSIoperator->uFESpace().mesh()->meanH() << std::endl
                  << "Solid mesh maxH     = " << M_FSIoperator->dFESpace().mesh()->maxH() << std::endl
                  << "Solid mesh meanH    = " << M_FSIoperator->dFESpace().mesh()->meanH() << std::endl << std::endl;
    }
}


// ===================================================
// MultiscaleInterfaceFluid Methods
// ===================================================
void
MultiscaleModelFSI3D::imposeBoundaryFlowRate( const bcFlag_Type& flag, const function_Type& function ) const
{
    BCFunctionBase base;
    base.setFunction( function );

    M_fluidBC->handler()->addBC( "CouplingFlowRate_Model_" + number2string( M_ID ) + "_Flag_" + number2string( flag ), flag, Flux, Full, base, 3 );
}

void
MultiscaleModelFSI3D::imposeBoundaryStress( const bcFlag_Type& flag, const function_Type& function ) const
{
    BCFunctionBase base;
    base.setFunction( function );

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
        return M_FSIoperator->fluid().lagrangeMultiplier(flag, *M_fluidBC->handler(), M_FSIoperator->solution() );
    else
        return M_FSIoperator->fluid().pressure( flag, M_FSIoperator->solution() );
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

    if ( multiscaleProblemStep > 0 )
    {
        M_importerFluid->setMeshProcId( M_FSIoperator->uFESpace().mesh(), M_FSIoperator->uFESpace().map().comm().MyPID() );
        M_importerSolid->setMeshProcId( M_FSIoperator->dFESpace().mesh(), M_FSIoperator->dFESpace().map().comm().MyPID() );

        M_importerFluid->addVariable( IOData_Type::VectorField, "Velocity (fluid)",     M_FSIoperator->uFESpacePtr(),  M_fluidVelocityAndPressure, UInt(0) );
        M_importerFluid->addVariable( IOData_Type::ScalarField, "Pressure (fluid)",     M_FSIoperator->pFESpacePtr(),  M_fluidVelocityAndPressure, UInt(3 * M_FSIoperator->uFESpace().dof().numTotalDof() ) );
        M_importerFluid->addVariable( IOData_Type::VectorField, "Displacement (fluid)", M_FSIoperator->mmFESpacePtr(), M_fluidDisplacement, UInt(0) );

        M_importerSolid->addVariable( IOData_Type::VectorField, "Velocity (solid)",     M_FSIoperator->dFESpacePtr(),  M_solidVelocity,     UInt(0) );
        M_importerSolid->addVariable( IOData_Type::VectorField, "Displacement (solid)", M_FSIoperator->dFESpacePtr(),  M_solidDisplacement, UInt(0) );

        // Import
        M_exporterFluid->setTimeIndex( M_importerFluid->importFromTime( M_data->dataFluid()->dataTime()->initialTime() ) + 1 );
        M_exporterSolid->setTimeIndex( M_importerSolid->importFromTime( M_data->dataSolid()->getdataTime()->initialTime() ) + 1 );

#ifdef HAVE_HDF5
        if ( M_FSIoperator->isFluid() )
            ( multiscaleDynamicCast< hdf5IOFile_Type >( M_importerFluid ) )->closeFile();
        if ( M_FSIoperator->isSolid() )
            ( multiscaleDynamicCast< hdf5IOFile_Type >( M_importerSolid ) )->closeFile();
#endif

        // Assemble the Monolithic solution
        vector_Type solution( *M_FSIoperator->couplingVariableMap() );

        // Add fluid
        vector_Type temporaryVector( *M_fluidVelocityAndPressure, Unique, Zero );
        solution = temporaryVector;

        // Add solid
        UInt offset = boost::dynamic_pointer_cast< FSIMonolithic > ( M_FSIoperator )->getOffset();

        temporaryVector = 0;
        temporaryVector.subset( *M_solidDisplacement, M_solidDisplacement->map(), static_cast<UInt> ( 0 ), offset );
        temporaryVector /= M_FSIoperator->solid().getRescaleFactor() * M_data->dataFluid()->dataTime()->timeStep();

        solution += temporaryVector;

        // Add harmonic extension
        if ( !M_data->method().compare("monolithicGI") )
        {
            offset = boost::dynamic_pointer_cast< FSIMonolithicGI > ( M_FSIoperator )->mapWithoutMesh().map( Unique )->NumGlobalElements();

            temporaryVector = 0;
            temporaryVector.subset( *M_fluidDisplacement, M_fluidDisplacement->map(), static_cast<UInt> ( 0 ), offset );

            solution += temporaryVector;
        }

        *M_fluidVelocityAndPressure = solution;
    }
    else
    {
        // Initialize solution
        *M_fluidVelocityAndPressure = 0.0;
        *M_fluidDisplacement        = 0.0;
        *M_solidVelocity            = 0.0;
        *M_solidDisplacement        = 0.0;

        // We initialize the fluid pressure equal to the external pressure
        vector_Type fluidPressure( M_FSIoperator->pFESpace().mapPtr(), Unique );
        fluidPressure = M_data->dataSolid()->externalPressure();

        vector_Type extendedFluidPressure( *M_fluidVelocityAndPressure, Unique, Zero );
        extendedFluidPressure.subset( fluidPressure, fluidPressure.map(), static_cast <UInt> ( 0 ), UInt( 3 * M_FSIoperator->uFESpace().dof().numTotalDof() ) );

        *M_fluidVelocityAndPressure += extendedFluidPressure;
    }

    // Initialize solution at time tn
    M_fluidVelocityAndPressure_tn.reset( new vector_Type( *M_fluidVelocityAndPressure ) );
    M_fluidDisplacement_tn.reset( new vector_Type( *M_fluidDisplacement ) );
    M_solidVelocity_tn.reset( new vector_Type( *M_solidVelocity ) );
    M_solidDisplacement_tn.reset( new vector_Type( *M_solidDisplacement ) );

    // Initialize all the quantities in the solver to time tn
    M_FSIoperator->initialize( M_fluidVelocityAndPressure_tn, M_fluidDisplacement_tn, M_solidVelocity_tn, M_solidDisplacement_tn );
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
MultiscaleModelFSI3D::updateSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::updateSolution() \n";
#endif

    if ( M_FSIoperator->isFluid() )
    {
        M_FSIoperator->exportFluidVelocityAndPressure( *M_fluidVelocityAndPressure );
        M_FSIoperator->exportFluidDisplacement( *M_fluidDisplacement );
    }

    if ( M_FSIoperator->isSolid() )
    {
        M_FSIoperator->exportSolidVelocity( *M_solidVelocity );
        M_FSIoperator->exportSolidDisplacement( *M_solidDisplacement );
    }
}

void
MultiscaleModelFSI3D::setupExporter( IOFilePtr_Type& exporter, const GetPot& dataFile, const std::string& label )
{
    const std::string exporterType = dataFile( "exporter/type", "ensight" );
#ifdef HAVE_HDF5
    if ( exporterType.compare( "hdf5" ) == 0 )
        exporter.reset( new hdf5IOFile_Type( dataFile, label ) );
    else
#endif
        exporter.reset( new ensightIOFile_Type( dataFile, label ) );

    exporter->setPrefix( "Step_" + number2string( multiscaleProblemStep ) + "_Model_" + number2string( M_ID ) + label );
    exporter->setPostDir( multiscaleProblemFolder );
}

void
MultiscaleModelFSI3D::setupImporter( IOFilePtr_Type& importer, const GetPot& dataFile, const std::string& label )
{
    const std::string importerType = dataFile( "exporter/type", "ensight" );
#ifdef HAVE_HDF5
    if ( importerType.compare( "hdf5" ) == 0 )
        importer.reset( new hdf5IOFile_Type( dataFile, label ) );
    else
#endif
        importer.reset( new ensightIOFile_Type( dataFile, label ) );

    importer->setPrefix( "Step_" + number2string( multiscaleProblemStep - 1 ) + "_Model_" + number2string( M_ID ) + label );
    importer->setPostDir( multiscaleProblemFolder );
}

void
MultiscaleModelFSI3D::setExporterFluid( const IOFilePtr_Type& exporter )
{
    M_fluidVelocityAndPressure.reset( new vector_Type( M_FSIoperator->fluid().getMap(),  M_exporterFluid->mapType() ) );
    M_fluidDisplacement.reset       ( new vector_Type( M_FSIoperator->mmFESpace().map(), M_exporterFluid->mapType() ) );

    exporter->setMeshProcId( M_FSIoperator->uFESpace().mesh(), M_FSIoperator->uFESpace().map().comm().MyPID() );

    exporter->addVariable( IOData_Type::VectorField, "Velocity (fluid)",     M_FSIoperator->uFESpacePtr(),  M_fluidVelocityAndPressure, static_cast<UInt> ( 0 ) );
    exporter->addVariable( IOData_Type::ScalarField, "Pressure (fluid)",     M_FSIoperator->pFESpacePtr(),  M_fluidVelocityAndPressure, static_cast<UInt> (3 * M_FSIoperator->uFESpace().dof().numTotalDof() ) );
    exporter->addVariable( IOData_Type::VectorField, "Displacement (fluid)", M_FSIoperator->mmFESpacePtr(), M_fluidDisplacement,        static_cast<UInt> ( 0 ) );
}

void
MultiscaleModelFSI3D::setExporterSolid( const IOFilePtr_Type& exporter )
{
    M_solidDisplacement.reset( new vector_Type( M_FSIoperator->dFESpace().map(), M_exporterSolid->mapType() ) );
    M_solidVelocity.reset    ( new vector_Type( M_FSIoperator->dFESpace().map(), M_exporterSolid->mapType() ) );

    exporter->setMeshProcId( M_FSIoperator->dFESpace().mesh(), M_FSIoperator->dFESpace().map().comm().MyPID() );

    exporter->addVariable( IOData_Type::VectorField, "Velocity (solid)",     M_FSIoperator->dFESpacePtr(), M_solidVelocity,     static_cast<UInt> ( 0 ) );
    exporter->addVariable( IOData_Type::VectorField, "Displacement (solid)", M_FSIoperator->dFESpacePtr(), M_solidDisplacement, static_cast<UInt> ( 0 ) );
}

void
MultiscaleModelFSI3D::setupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::setupLinearModel() \n";
#endif

    // The linear BCHandler is a copy of the original BCHandler with all BCFunctions giving zero
    bcPtr_Type linearBCHandler ( new bc_Type( *M_fluidBC->handler() ) );
    M_linearBC = linearBCHandler;

    // Set all the BCFunctions to zero
    BCFunctionBase bcBaseDeltaZero;
    bcBaseDeltaZero.setFunction( boost::bind( &MultiscaleModelFSI3D::bcFunctionDeltaZero, this, _1, _2, _3, _4, _5 ) );

    for ( bc_Type::bcBaseIterator_Type i = M_linearBC->begin() ; i != M_linearBC->end() ; ++i )
        i->setBCFunction( bcBaseDeltaZero );

    // Setup linear solution & the RHS
    M_linearSolution.reset( new vector_Type( M_FSIoperator->un()->map() ) );
    M_linearRHS.reset( new vector_Type( M_FSIoperator->un()->map() ) );
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
