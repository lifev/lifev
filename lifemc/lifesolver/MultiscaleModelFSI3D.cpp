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
 *  @brief File containing the MultiScale Model FSI3D
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
        M_FSIoperator                  (),
        M_data                         ( new data_Type() ),
        M_exporterFluid                (),
        M_exporterSolid                (),
        M_importerFluid                (),
        M_importerSolid                (),
        M_fluidVelocityPressure        (),
        M_fluidDisplacement            (),
        M_solidVelocity                (),
        M_solidDisplacement            (),
        M_fluidVelocityPressure_tn     (),
        M_solidDisplacement_tn         (),
        M_solidDisplacementOld_tn      (),
        M_rhs_tn                       (),
        M_nonLinearRichardsonIteration (),
        M_fluidBC                      ( new bcInterface_Type() ),
        M_solidBC                      ( new bcInterface_Type() ),
        M_harmonicExtensionBC          ( new bcInterface_Type() ),
        M_linearBC                     (),
        M_linearRHS                    (),
        M_linearSolution               (),
        M_bcBaseDeltaZero              (),
        M_bcBaseDeltaOne               ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::MultiscaleModelFSI3D() \n";
#endif

    M_type = FSI3D;

    BlockPrecFactory::instance().registerProduct("ComposedDNND",      &ComposedDNND::createComposedDNND);
    BlockPrecFactory::instance().registerProduct("AdditiveSchwarz",   &BlockMatrix::createAdditiveSchwarz) ;
    BlockPrecFactory::instance().registerProduct("AdditiveSchwarzRN", &BlockMatrixRN::createAdditiveSchwarzRN ) ;
    BlockPrecFactory::instance().registerProduct("ComposedDN",        &ComposedDN::createComposedDN ) ;
    BlockPrecFactory::instance().registerProduct("ComposedDN2",       &ComposedDN::createComposedDN2 );

    BlockMatrix::Factory::instance().registerProduct("AdditiveSchwarz",   &BlockMatrix::createAdditiveSchwarz ) ;
    BlockMatrix::Factory::instance().registerProduct("AdditiveSchwarzRN", &BlockMatrixRN::createAdditiveSchwarzRN ) ;

    FSIOperator_Type::solid_Type::StructureSolverFactory::instance().registerProduct( "linearVenantKirchhof", &FSIOperator_Type::createLinearStructure );
//    FSIOperator_Type::solid_Type::StructureSolverFactory::instance().registerProduct( "nonLinearVenantKirchhof", &FSIOperator_Type::createNonLinearStructure );
}

// ===================================================
// MultiScale PhysicalModel Virtual Methods
// ===================================================
void
MultiscaleModelFSI3D::setupData( const std::string& fileName )
{
    multiscaleModel_Type::setupData( fileName );

    GetPot dataFile( fileName );

    // Load data
    M_data->setup( dataFile );

    if ( M_globalData.get() )
        setupGlobalData( fileName );

    // Create FSIOperator
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
    Int numLM = M_FSIoperator->imposeFlux();
    M_FSIoperator->setFluxesNumber( numLM );
    M_FSIoperator->setupFluidSolid( numLM );

    // Setup Exporters
    if ( M_FSIoperator->isFluid() )
        setExporterFluid( M_exporterFluid );

    if ( M_FSIoperator->isSolid() )
        setExporterSolid( M_exporterSolid );

    //Setup solution
    initializeSolution();

    //Setup linear model
    setupLinearModel();

    M_solidDisplacement_tn.reset( new vector_Type( *M_fluidDisplacement ));
}

void
MultiscaleModelFSI3D::buildSystem()
{
    // Update BCInterface Operator BC
    updateBC();

    M_FSIoperator->setupSystem();
    M_FSIoperator->buildSystem();
    M_FSIoperator->updateSystem();

    // Update solution at time n
    M_fluidVelocityPressure_tn.reset( new vector_Type( M_FSIoperator->solution() ) );
    M_rhs_tn.reset( new vector_Type( M_FSIoperator->solution().map() ) );
    M_solidDisplacement_tn.reset( new vector_Type( M_FSIoperator->meshDisp() ) );
    M_solidDisplacementOld_tn.reset( new vector_Type( M_FSIoperator->meshDisp() ) );
}

void
MultiscaleModelFSI3D::updateSystem()
{
    // Update BCInterface Operator BC
    updateBC();

    // Update System
    M_FSIoperator->updateSystem();

    //Update solution at time n
    *M_fluidVelocityPressure_tn = *M_FSIoperator->solutionPtr();
    *M_solidDisplacementOld_tn = M_FSIoperator->meshMotion().dispOld();
    *M_solidDisplacement_tn = M_FSIoperator->meshDisp();
    *M_rhs_tn *= 0;

    M_nonLinearRichardsonIteration = 0;

    boost::dynamic_pointer_cast<Monolithic>(M_FSIoperator)->precPtrView()->setRecompute( 1, true );
}

void
MultiscaleModelFSI3D::solveSystem( )
{
    UInt maxSubIterationNumber = M_data->maxSubIterationNumber();
    std::ofstream outRes; // Unuseful variable

    // Non-linear Richardson solver
    vectorPtr_Type solution( new vector_Type( *M_fluidVelocityPressure_tn ) ) ;
    boost::dynamic_pointer_cast<LifeV::Monolithic>(M_FSIoperator)->initializeMesh(M_solidDisplacement_tn);
    M_FSIoperator->fluid().initialize( *M_fluidVelocityPressure );//useless?
    vectorPtr_Type newDisp(new vector_Type(*M_solidDisplacement));
    M_FSIoperator->solid().initialize( newDisp, M_solidVelocity );
    vectorPtr_Type newRHS(new vector_Type(*M_rhs_tn));
    M_FSIoperator->setRHS( newRHS );
    M_FSIoperator->setSolution( *M_fluidVelocityPressure_tn );
    M_FSIoperator->initializeBDF( *M_fluidVelocityPressure_tn );
    M_FSIoperator->setRestarts( true );

    if (M_nonLinearRichardsonIteration != 0)
    {
        boost::dynamic_pointer_cast<Monolithic>(M_FSIoperator)->precPtrView()->setRecompute( 1, false );
        M_FSIoperator->updateRHS();
        M_FSIoperator->applyBoundaryConditions( );
    }

    UInt status = nonLinRichardson( *solution, *M_FSIoperator,
                                    M_data->absoluteTolerance(),
                                    M_data->relativeTolerance(),
                                    maxSubIterationNumber,
                                    M_data->errorTolerance(),
                                    M_data->linesearch(),
                                    outRes,
                                    M_data->dataFluid()->dataTime()->time(),
                                    M_nonLinearRichardsonIteration
                                  );
    if (M_nonLinearRichardsonIteration == 0)
        *M_fluidDisplacement = M_FSIoperator->meshDisp();

    M_nonLinearRichardsonIteration=1;

    // We set the solution for the next time step
    M_FSIoperator->setSolutionPtr( solution );

    if ( status == EXIT_FAILURE )
        std::cout << "Non-Linear Richardson failed to converge" << std::endl;
}

void
MultiscaleModelFSI3D::saveSolution()
{
    // TODO Post-process must be made here. We need to add to HDF5exporter some methods to:
    // 1) save only the xmf file (done once, as it is independent from the variable)
    // 2) save a specific variable at a specific time (defined by variablename && time)
    /*
        // save xmf file for the fluid at the  present time step
        if( M_FSIoperator->isFluid() )
        {
            // if ( segregated && semiimplicit false ) || fullMonolithic
            // save fluidDisplacement at this time step

            // if ( segregated && semiimplicit true)  || monolithic
            // save fluidDisplacement at the previous time step

            // save velocity and pressure always at this time step
        }

        // For the solid we can use classical approach
        if( M_FSIoperator->isSolid() )
        {
            M_FSIoperator->getSolidDisp( *M_solidDisplacement );// displacement(), M_offset);
            M_FSIoperator->getSolidVel( *M_solidVelocity );// displacement(), M_offset);
        }
    */

    // TODO Post-processing is not working with MonolithicGI + it is also not saving last time step for MonolithicGE
    //      To solve this do HE after FluidAndSolid

    // Saving Fluid (displacement) for this post-processing
//     if( M_FSIoperator->isFluid() )
//         *M_fluidDisplacement = M_FSIoperator->meshDisp();

    if ( M_FSIoperator->isFluid() )
        M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->time() - M_data->dataFluid()->dataTime()->timeStep() );
    if ( M_FSIoperator->isSolid() )
        M_exporterSolid->postProcess( M_data->dataSolid()->dataTime()->time() - M_data->dataSolid()->dataTime()->timeStep() );

#ifdef HAVE_HDF5
    if ( M_data->dataFluid()->dataTime()->isLastTimeStep() )
    {
        if ( M_FSIoperator->isFluid() )
            ( multiscaleDynamicCast< hdf5IOFile_Type >( M_exporterFluid ) )->closeFile();
        if ( M_FSIoperator->isSolid() )
            ( multiscaleDynamicCast< hdf5IOFile_Type >( M_exporterSolid ) )->closeFile();
    }
#endif

    // Saving Fluid (velocity and pressure) and Solid (displacement) for next post-processing
    if ( M_FSIoperator->isFluid() )
        M_FSIoperator->getFluidVelAndPres( *M_fluidVelocityPressure );
    if ( M_FSIoperator->isSolid() )
    {
        M_FSIoperator->getSolidDisp( *M_solidDisplacement );
        M_FSIoperator->getSolidVel( *M_solidVelocity );
    }
}

void
MultiscaleModelFSI3D::showMe()
{
    if ( M_displayer->isLeader() )
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

        std::cout << "Fluid mesh maxH     = " << M_FSIoperator->uFESpace().mesh()->maxH() << std::endl
                  << "Fluid mesh meanH    = " << M_FSIoperator->uFESpace().mesh()->meanH() << std::endl
                  << "Solid mesh maxH     = " << M_FSIoperator->dFESpace().mesh()->maxH() << std::endl
                  << "Solid mesh meanH    = " << M_FSIoperator->dFESpace().mesh()->meanH() << std::endl << std::endl;
    }
}

// ===================================================
// Methods
// ===================================================
void
MultiscaleModelFSI3D::setupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MS_Model_Fluid3D::SetupLinearModel( ) \n";
#endif

    // Define BCFunctions for tangent problem
    M_bcBaseDeltaZero.setFunction( boost::bind( &MultiscaleModelFSI3D::bcFunctionDeltaZero, this, _1, _2, _3, _4, _5 ) );
    M_bcBaseDeltaOne.setFunction(  boost::bind( &MultiscaleModelFSI3D::bcFunctionDeltaOne,  this, _1, _2, _3, _4, _5 ) );

    // The linear BCHandler is a copy of the original BCHandler with all BCFunctions giving zero
    bcPtr_Type LinearBCHandler ( new bc_Type( *M_fluidBC->handler() ) );
    M_linearBC = LinearBCHandler;

    // Set all the BCFunctions to zero
    for ( bc_Type::bcBaseIterator_Type i = M_linearBC->begin() ; i != M_linearBC->end() ; ++i )
        i->setBCFunction( M_bcBaseDeltaZero );

    // Setup linear solution & the RHS
    M_linearSolution.reset( new vector_Type( M_FSIoperator->un()->map() ) );
    M_linearRHS.reset( new vector_Type( M_FSIoperator->un()->map() ) );
}

void
MultiscaleModelFSI3D::updateLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::UpdateLinearModel() \n";
#endif

    //Create the RHS
    *M_linearRHS *= 0;
    M_FSIoperator->bcManageVectorRHS( M_linearBC, *M_linearRHS );
}

void
MultiscaleModelFSI3D::solveLinearModel( bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::SolveLinearModel() \n";
#endif

    if ( !solveLinearSystem )
        return;

    imposePerturbation();

    updateLinearModel();

    //Solve the linear problem
    M_FSIoperator->solveJac( *M_linearSolution, *M_linearRHS, 0. );

    resetPerturbation();

    //This flag avoid recomputation of the same system
    solveLinearSystem = false;
}

// ===================================================
// Get Methods
// ===================================================
Real
MultiscaleModelFSI3D::boundaryStress( const bcFlag_Type& flag, const stress_Type& stressType ) const
{
    switch ( stressType )
    {
    case StaticPressure:
    {
        return -boundaryPressure( flag );
    }

    case TotalPressure:
    {
        return -boundaryPressure( flag ) + boundaryDynamicPressure( flag ) * ( ( boundaryFlowRate( flag ) > 0.0 ) ? 1 : -1 );
    }

    case LagrangeMultiplier:
    {
        return -boundaryLagrangeMultiplier( flag );
    }

    default:

        std::cout << "ERROR: Invalid stress type [" << Enum2String( stressType, multiscaleStressesMap ) << "]" << std::endl;

        return 0.0;
    }
}

Real
MultiscaleModelFSI3D::boundaryDeltaFlowRate( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    solveLinearModel( solveLinearSystem );

    return M_FSIoperator->fluid().flux( flag, *M_linearSolution );
}

Real
MultiscaleModelFSI3D::boundaryDeltaPressure( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    solveLinearModel( solveLinearSystem );

    return M_FSIoperator->fluid().pressure( flag, *M_linearSolution );
}

Real
MultiscaleModelFSI3D::boundaryDeltaDynamicPressure( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    return boundaryDensity( flag ) * boundaryDeltaFlowRate( flag, solveLinearSystem ) * boundaryFlowRate( flag ) / ( boundaryArea( flag ) * boundaryArea( flag ) );
}

Real
MultiscaleModelFSI3D::boundaryDeltaLagrangeMultiplier( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    solveLinearModel( solveLinearSystem );

    return M_FSIoperator->fluid().LagrangeMultiplier( flag, *M_linearBC, *M_linearSolution );
}

Real
MultiscaleModelFSI3D::boundaryDeltaStress( const bcFlag_Type& flag, bool& solveLinearSystem, const stress_Type& stressType )
{
    switch ( stressType )
    {
    case StaticPressure:
    {
        return -boundaryDeltaPressure( flag, solveLinearSystem );
    }

    case TotalPressure:
    {
        return -boundaryDeltaPressure( flag, solveLinearSystem ) + boundaryDeltaDynamicPressure( flag, solveLinearSystem ); //Verify the sign of DynamicPressure contribute!
    }

    case LagrangeMultiplier:
    {
        return -boundaryDeltaLagrangeMultiplier( flag, solveLinearSystem );
    }

    default:

        std::cout << "ERROR: Invalid stress type [" << Enum2String( stressType, multiscaleStressesMap ) << "]" << std::endl;

        return 0.0;
    }
}

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleModelFSI3D::setupGlobalData( const std::string& fileName )
{
    GetPot dataFile( fileName );

    //Global data time
    M_data->dataFluid()->setDataTime( M_globalData->dataTime() );
    M_data->dataSolid()->setDataTime( M_globalData->dataTime() );

    //Global physical quantities
    if ( !dataFile.checkVariable( "fluid/physics/density" ) )
        M_data->dataFluid()->density( M_globalData->fluidDensity() );
    if ( !dataFile.checkVariable( "fluid/physics/viscosity" ) )
        M_data->dataFluid()->viscosity( M_globalData->fluidViscosity() );

    if ( !dataFile.checkVariable( "solid/physics/density" ) )
        M_data->dataSolid()->setDensity( M_globalData->structureDensity() );
    if ( !dataFile.checkVariable( "solid/physics/poisson" ) )
        M_data->dataSolid()->setPoisson( M_globalData->structurePoissonCoefficient() );
    if ( !dataFile.checkVariable( "solid/physics/young" ) )
        M_data->dataSolid()->setYoung( M_globalData->structureYoungModulus() );

    //M_data->showMe();
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
    else
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
    M_fluidVelocityPressure.reset( new vector_Type( M_FSIoperator->fluid().getMap(),  M_exporterFluid->mapType() ) );
    M_fluidDisplacement.reset    ( new vector_Type( M_FSIoperator->mmFESpace().map(), M_exporterFluid->mapType() ) );

    exporter->setMeshProcId( M_FSIoperator->uFESpace().mesh(),
                             M_FSIoperator->uFESpace().map().comm().MyPID() );

    exporter->addVariable( ExporterData::Vector, "Fluid Velocity", M_fluidVelocityPressure, UInt(0), M_FSIoperator->uFESpace().dof().numTotalDof()  );
    exporter->addVariable( ExporterData::Scalar, "Fluid Pressure", M_fluidVelocityPressure, UInt(3 * M_FSIoperator->uFESpace().dof().numTotalDof()  ),
                           UInt(    M_FSIoperator->pFESpace().dof().numTotalDof()  ) );
    exporter->addVariable( ExporterData::Vector, "Fluid Displacement", M_fluidDisplacement, UInt(0), M_FSIoperator->mmFESpace().dof().numTotalDof() );
}

void
MultiscaleModelFSI3D::setExporterSolid( const IOFilePtr_Type& exporter )
{
    M_solidDisplacement.reset( new vector_Type( M_FSIoperator->dFESpace().map(), M_exporterSolid->mapType() ) );
    M_solidVelocity.reset    ( new vector_Type( M_FSIoperator->dFESpace().map(), M_exporterSolid->mapType() ) );

    exporter->setMeshProcId( M_FSIoperator->dFESpace().mesh(),
                             M_FSIoperator->dFESpace().map().comm().MyPID() );

    exporter->addVariable( ExporterData::Vector, "Solid Velocity",     M_solidVelocity,     UInt(0), M_FSIoperator->dFESpace().dof().numTotalDof() );
    exporter->addVariable( ExporterData::Vector, "Solid Displacement", M_solidDisplacement, UInt(0), M_FSIoperator->dFESpace().dof().numTotalDof() );
}

void
MultiscaleModelFSI3D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::InitializeSolution() \n";
#endif

    if ( multiscaleProblemStep > 0 )
    {
        M_importerFluid->setMeshProcId( M_FSIoperator->uFESpace().mesh(), M_FSIoperator->uFESpace().map().comm().MyPID() );
        M_importerSolid->setMeshProcId( M_FSIoperator->dFESpace().mesh(), M_FSIoperator->dFESpace().map().comm().MyPID() );

        M_importerFluid->addVariable( ExporterData::Vector, "Fluid Velocity", M_fluidVelocityPressure, UInt(0), M_FSIoperator->uFESpace().dof().numTotalDof()  );
        M_importerFluid->addVariable( ExporterData::Scalar, "Fluid Pressure", M_fluidVelocityPressure, UInt(3 * M_FSIoperator->uFESpace().dof().numTotalDof()  ),
                                      UInt(    M_FSIoperator->pFESpace().dof().numTotalDof()  ) );
        M_importerFluid->addVariable( ExporterData::Vector, "Fluid Displacement", M_fluidDisplacement, UInt(0), M_FSIoperator->mmFESpace().dof().numTotalDof() );

        M_importerSolid->addVariable( ExporterData::Vector, "Solid Velocity",     M_solidVelocity,     UInt(0), M_FSIoperator->dFESpace().dof().numTotalDof() );
        M_importerSolid->addVariable( ExporterData::Vector, "Solid Displacement", M_solidDisplacement, UInt(0), M_FSIoperator->dFESpace().dof().numTotalDof() );

        // Import
        M_exporterFluid->setStartIndex( M_importerFluid->importFromTime( M_data->dataFluid()->dataTime()->initialTime() ) + 1 );
        M_exporterSolid->setStartIndex( M_importerSolid->importFromTime( M_data->dataSolid()->dataTime()->initialTime() ) + 1 );

        // Set Initial solution
        // IMPORTANT NOTE:
        // 1) TODO Remove nRestart flag from the file
        // 2) TODO This part should be rewritten better
        vectorPtr_Type UniqueVFDOld( new vector_Type( *M_fluidDisplacement, Unique, Zero ) );
        dynamic_cast< Monolithic* >( M_FSIoperator.get())->initializeMesh( UniqueVFDOld );

        vectorPtr_Type initSol( new EpetraVector( *M_FSIoperator->couplingVariableMap() ) );
        vectorPtr_Type UniqueV( new vector_Type( *M_fluidVelocityPressure, Unique, Zero ) );
        *initSol = *UniqueV;
        M_FSIoperator->fluid().initialize( *initSol );

        UniqueV.reset( new vector_Type( *M_FSIoperator->couplingVariableMap(), Unique, Zero ) );
        UInt offset = dynamic_cast<Monolithic*>(M_FSIoperator.get())->getOffset();
        UniqueV->subset( *M_solidDisplacement, M_solidDisplacement->map(), (UInt)0, offset );
        *UniqueV *= 1 / ( M_FSIoperator->solid().getRescaleFactor() * M_data->dataFluid()->dataTime()->timeStep() );
        M_FSIoperator->solid().initialize( UniqueV );
        *initSol += *UniqueV;

        if ( !M_data->method().compare("monolithicGI") )
        {
            vectorPtr_Type UniqueVFD ( new vector_Type( *M_FSIoperator->couplingVariableMap(), Unique, Zero ) );
            UniqueVFD->subset( *M_fluidDisplacement, M_fluidDisplacement->map(), (UInt)0, dynamic_cast<MonolithicGI*>(M_FSIoperator.get())->mapWithoutMesh().map(Unique)->NumGlobalElements());
            *initSol += *UniqueVFD;
        }

        vectorPtr_Type initSolSVel( new vector_Type( *M_FSIoperator->couplingVariableMap(), Unique, Zero ) );
        initSolSVel->subset( *M_solidVelocity,M_solidVelocity->map(), (UInt)0, offset );
        *initSolSVel *= 1 / ( M_FSIoperator->solid().getRescaleFactor() * M_data->dataSolid()->dataTime()->timeStep() );
        M_FSIoperator->solid().initializeVel( *initSolSVel );

        M_FSIoperator->setSolution( *initSol );
        M_FSIoperator->initializeBDF( *initSol );
    }
    else
    {
        M_FSIoperator->setSolution( *M_fluidVelocityPressure );
        M_FSIoperator->initializeBDF( *M_fluidVelocityPressure );
    }
}

void
MultiscaleModelFSI3D::imposePerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MS_Model_FSID::ImposePerturbation() \n";
#endif

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            M_linearBC->findBCWithFlag( ( *i )->flag( ( *i )->modelGlobalToLocalID( M_ID ) ) ).setBCFunction( M_bcBaseDeltaOne );

            break;
        }
}

void
MultiscaleModelFSI3D::resetPerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MultiscaleModelFSI3D::ResetPerturbation() \n";
#endif

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            M_linearBC->findBCWithFlag( ( *i )->flag( ( *i )->modelGlobalToLocalID( M_ID ) ) ).setBCFunction( M_bcBaseDeltaZero );

            break;
        }
}

} // Namespace multiscale
} // Namespace LifeV
