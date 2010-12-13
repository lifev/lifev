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
 *  @brief File containing the MultiScale Model Fluid3D
 *
 *  @date 12-03-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleModelFluid3D.hpp>

namespace LifeV
{
namespace multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModelFluid3D::MultiscaleModelFluid3D() :
        multiscaleModel_Type           (),
        M_exporter                     (),
        M_importer                     (),
        M_fileName                     (),
        M_fluid                        (),
        M_bc                           ( new bcInterface_Type() ),
        M_bdf                          (),
        M_data                         ( new data_Type() ),
        M_dataMesh                     ( new DataMesh()),
        M_mesh                         (),
        M_map                          (),
        M_solution                     (),
        M_linearBC                     ( new bc_Type() ),
        M_updateLinearModel            ( true ),
        M_uFESpace                     (),
        M_pFESpace                     (),
        M_lmDOF                        ( 0 ),
        M_alpha                        ( 0 ),
        M_beta                         (),
        M_rhs                          (),
        M_subiterationsMaximumNumber   (),
        M_tolerance                    (),
        M_generalizedAitken            ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::MultiscaleModelFluid3D() \n";
#endif

    M_type = Fluid3D;
}

// ===================================================
// MultiScale PhysicalModel Virtual Methods
// ===================================================
void
MultiscaleModelFluid3D::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::SetupData( ) \n";
#endif

    multiscaleModel_Type::setupData( fileName );
    M_fileName = fileName;

    GetPot dataFile( fileName );

    //Fluid data
    M_data->setup( dataFile );
    if ( M_globalData.get() )
        setupGlobalData( fileName );

    M_dataMesh->setup(dataFile, "fluid/space_discretization");

    // Parameters for the NS Iterations
    M_subiterationsMaximumNumber = dataFile( "fluid/miscellaneous/SubITMax", 0 );
    M_tolerance                  = dataFile( "fluid/miscellaneous/Tolerance", 1.e-6 );

    M_generalizedAitken.setDefaultOmega(     dataFile( "fluid/miscellaneous/Omega",        1.e-3 ) );
    M_generalizedAitken.setOmegaMin(         dataFile( "fluid/miscellaneous/range",        M_generalizedAitken.defaultOmegaFluid()/1024, 0 ) );
    M_generalizedAitken.setOmegaMax(         dataFile( "fluid/miscellaneous/range",        M_generalizedAitken.defaultOmegaFluid()*1024, 1 ) );
    M_generalizedAitken.useDefaultOmega(     dataFile( "fluid/miscellaneous/fixedOmega",   false ) );
    M_generalizedAitken.setMinimizationType( dataFile( "fluid/miscellaneous/inverseOmega", true ) );

    //Boundary Conditions for the problem
    M_bc->fillHandler( fileName, "fluid" );

    //Setup Exporter & Importer
    setupExporterImporter( fileName );
}

void
MultiscaleModelFluid3D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::SetupProblem() \n";
#endif

    //Mesh
    setupMesh();

    //FEspace
    setupFEspace();

    //Add flow rate offset to BC
    M_lmDOF = M_bc->handler()->getNumberBCWithType( Flux );
    setupBCOffset( M_bc->handler() );

    //Fluid
    M_fluid.reset( new fluid_Type( M_data, *M_uFESpace, *M_pFESpace, M_comm, M_lmDOF ) );
    M_bc->setPhysicalSolver( M_fluid );

    GetPot dataFile( M_fileName );
    M_fluid->setUp( dataFile ); //Remove Preconditioner and Solver if possible!

    //Fluid MAP
    M_map.reset( new EpetraMap( M_fluid->getMap() ) );

    //BDF
    M_bdf.reset( new bdf_Type( M_data->dataTime()->getBDF_order() ) );

    //Problem coefficients
    M_beta.reset( new fluidVector_Type( M_map ) );
    M_rhs.reset ( new fluidVector_Type( M_map ) );

    //Post-processing
    M_exporter->setMeshProcId( M_mesh->mesh(), M_comm->MyPID() );

    M_solution.reset( new fluidVector_Type( *M_fluid->solution(), M_exporter->mapType() ) );
    if ( M_exporter->mapType() == Unique )
        M_solution->setCombineMode( Zero );

    M_exporter->addVariable( ExporterData::Vector, "Fluid Velocity", M_solution, static_cast<UInt> ( 0 ), M_uFESpace->dof().numTotalDof() );
    M_exporter->addVariable( ExporterData::Scalar, "Fluid Pressure", M_solution, 3 * M_uFESpace->dof().numTotalDof(),              M_pFESpace->dof().numTotalDof() );

    //Setup linear model
    setupLinearModel();

    //Setup solution
    initializeSolution();
}

void
MultiscaleModelFluid3D::buildSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::BuildSystem() \n";
#endif

    //Build constant matrices
    M_fluid->buildSystem();

    //Initialize BDF
    M_bdf->bdf_u().initialize_unk( *M_fluid->solution() );

    //Define problem coefficients
    if ( M_data->Stokes() )
    {
        M_alpha  = 0.0;
        *M_beta  = *M_fluid->solution(); //It is a stationary Navier-Stokes
        *M_rhs  *= 0.0;
    }
    else
    {
        M_alpha = M_bdf->bdf_u().coeff_der( 0 ) / M_data->dataTime()->getTimeStep();
        *M_beta = M_bdf->bdf_u().extrap();
        *M_rhs  = M_fluid->matrMass() * M_bdf->bdf_u().time_der( M_data->dataTime()->getTimeStep() );
    }

    //Set problem coefficients
    M_fluid->updateSystem( M_alpha, *M_beta, *M_rhs );
}

void
MultiscaleModelFluid3D::updateSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::UpdateSystem() \n";
#endif

    //Update BDF
    M_bdf->bdf_u().shift_right( *M_fluid->solution() );

    //Update problem coefficients
    M_alpha = M_bdf->bdf_u().coeff_der( 0 ) / M_data->dataTime()->getTimeStep();
    *M_beta = M_bdf->bdf_u().extrap();
    *M_rhs  = M_fluid->matrMass() * M_bdf->bdf_u().time_der( M_data->dataTime()->getTimeStep() );

    //Set problem coefficients
    M_fluid->updateSystem( M_alpha, *M_beta, *M_rhs );

    //Update operator BC
    M_bc->updatePhysicalSolverVariables();

    //Recompute preconditioner
    M_fluid->resetPrec( true );

    //Linear system need to be updated
    M_updateLinearModel = true;
}

void
MultiscaleModelFluid3D::solveSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::SolveSystem() \n";
#endif

    //Solve the problem
    M_fluid->iterate( *M_bc->handler() );

    if ( M_subiterationsMaximumNumber > 0 )
    {
        Real residual = ( *M_beta - *M_fluid->solution() ).Norm2(); // residual is computed on the whole solution vector;

        if ( M_displayer->isLeader() )
            std::cout << "  F-  Residual:                                " << residual << std::endl;

        M_generalizedAitken.restart();
        for ( UInt subIT = 1; subIT <= M_subiterationsMaximumNumber; ++subIT )
        {
            // Verify tolerance
            if ( residual <= M_tolerance )
                break;

            *M_beta += M_generalizedAitken.computeDeltaLambdaScalar( *M_beta, *M_beta - *M_fluid->solution() );

            //Linear model need to be updated!
            M_fluid->updateSystem( M_alpha, *M_beta, *M_rhs );
            M_bc->updatePhysicalSolverVariables();
            M_updateLinearModel = true;

            //Solve system
            M_fluid->iterate( *M_bc->handler() );

            // Check the new residual
            residual = ( *M_beta - *M_fluid->solution() ).Norm2(); // residual is computed on the whole solution vector

            // Display subiteration information
            if ( M_displayer->isLeader() )
            {
                std::cout << "  F-  Sub-iteration n.:                        " << subIT << std::endl;
                std::cout << "  F-  Residual:                                " << residual << std::endl;
            }
        }
    }
}

void
MultiscaleModelFluid3D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::SaveSolution() \n";
#endif

    //Post-processing
    *M_solution = *M_fluid->solution();
    M_exporter->postProcess( M_data->dataTime()->getTime() );

#ifdef HAVE_HDF5
    if ( M_data->dataTime()->isLastTimeStep() )
        ( multiscaleDynamicCast< hdf5IOFile_Type >( M_exporter ) )->CloseFile();
#endif

}

void
MultiscaleModelFluid3D::showMe()
{
    if ( M_displayer->isLeader() )
    {
        multiscaleModel_Type::showMe();

        std::cout << "Velocity FE order   = " << M_data->uOrder() << std::endl
                  << "Pressure FE order   = " << M_data->pOrder() << std::endl << std::endl;

        std::cout << "Velocity DOF        = " << 3 * M_uFESpace->dof().numTotalDof() << std::endl
                  << "Pressure DOF        = " << M_pFESpace->dof().numTotalDof() << std::endl
                  << "lmDOF               = " << M_lmDOF << std::endl << std::endl;

        std::cout << "Fluid mesh maxH     = " << M_mesh->mesh()->maxH() << std::endl
                  << "Fluid mesh meanH    = " << M_mesh->mesh()->meanH() << std::endl << std::endl;

        std::cout << "NS SubITMax         = " << M_subiterationsMaximumNumber << std::endl
                  << "NS Tolerance        = " << M_tolerance << std::endl << std::endl << std::endl << std::endl;
    }
}


// ===================================================
// Methods
// ===================================================
void
MultiscaleModelFluid3D::setupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::SetupLinearModel( ) \n";
#endif

    // Define BCFunctions for tangent problem
    M_bcBaseDeltaZero.setFunction( boost::bind( &MultiscaleModelFluid3D::bcFunctionDeltaZero, this, _1, _2, _3, _4, _5 ) );
    M_bcBaseDeltaOne.setFunction(  boost::bind( &MultiscaleModelFluid3D::bcFunctionDeltaOne,  this, _1, _2, _3, _4, _5 ) );

    // The linear BCHandler is a copy of the original BCHandler with all BCFunctions giving zero
    bcPtr_Type LinearBCHandler ( new bc_Type( *M_bc->handler() ) );
    M_linearBC = LinearBCHandler;

    // Set all te BCFunctions to zero
    for ( bc_Type::BCBase_Iterator i = M_linearBC->begin() ; i != M_linearBC->end() ; ++i )
        i->setBCFunction( M_bcBaseDeltaZero );
}

void
MultiscaleModelFluid3D::updateLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::UpdateLinearModel() \n";
#endif

    //Create an empty vector
    fluidVector_Type vectorZero( *M_solution );
    vectorZero = 0.0;

    //updateLinearModel TODO REMOVE ?
    M_fluid->updateLinearSystem( M_fluid->matrNoBC(),
                                 M_alpha,
                                 *M_beta,
                                 *M_fluid->solution(),
                                 vectorZero,
                                 vectorZero,
                                 vectorZero,
                                 vectorZero );

    //Linear System Updated
    M_updateLinearModel = false;
}

void
MultiscaleModelFluid3D::solveLinearModel( bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::SolveLinearModel() \n";
#endif

    if ( !solveLinearSystem )
        return;

    imposePerturbation();

    if ( M_updateLinearModel )
        updateLinearModel();

    //Solve the linear problem
    M_fluid->iterateLin( *M_linearBC );

    resetPerturbation();

    //This flag avoid recomputation of the same system
    solveLinearSystem = false;
}

// ===================================================
// Set Methods
// ===================================================
void
MultiscaleModelFluid3D::setSolution( const fluidVectorPtr_Type& solution )
{
    M_solution = solution;

    M_fluid->initialize( *M_solution );
}

// ===================================================
// Get Methods (couplings)
// ===================================================
Real
MultiscaleModelFluid3D::boundaryStress( const BCFlag& flag, const stress_Type& stressType ) const
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
MultiscaleModelFluid3D::boundaryDeltaFlowRate( const BCFlag& flag, bool& solveLinearSystem )
{
    solveLinearModel( solveLinearSystem );

    return M_fluid->GetLinearFlux( flag );
}

Real
MultiscaleModelFluid3D::boundaryDeltaPressure( const BCFlag& flag, bool& solveLinearSystem )
{
    solveLinearModel( solveLinearSystem );

    return M_fluid->GetLinearPressure( flag );
}

Real
MultiscaleModelFluid3D::boundaryDeltaDynamicPressure( const BCFlag& flag, bool& solveLinearSystem )
{
    return boundaryDensity( flag ) * boundaryDeltaFlowRate( flag, solveLinearSystem ) * boundaryFlowRate( flag ) / ( boundaryArea( flag ) * boundaryArea( flag ) );
}

Real
MultiscaleModelFluid3D::boundaryDeltaLagrangeMultiplier( const BCFlag& flag, bool& solveLinearSystem )
{
    solveLinearModel( solveLinearSystem );

    return M_fluid->LinearLagrangeMultiplier( flag, *M_linearBC );
}

Real
MultiscaleModelFluid3D::boundaryDeltaStress( const BCFlag& flag, bool& solveLinearSystem, const stress_Type& stressType )
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
MultiscaleModelFluid3D::setupGlobalData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::SetupGlobalData( fileName ) \n";
#endif

    GetPot dataFile( fileName );

    //Global data time
    M_data->setDataTime( M_globalData->dataTime() );

    //Global physical quantities
    if ( !dataFile.checkVariable( "fluid/physics/density" ) )
        M_data->density( M_globalData->fluidDensity() );
    if ( !dataFile.checkVariable( "fluid/physics/viscosity" ) )
        M_data->viscosity( M_globalData->fluidViscosity() );
}

void
MultiscaleModelFluid3D::setupExporterImporter( const std::string& fileName )
{
    GetPot dataFile( fileName );

    //Exporter
    const std::string exporterType = dataFile( "exporter/type", "ensight" );

#ifdef HAVE_HDF5
    if ( !exporterType.compare( "hdf5" ) )
        M_exporter.reset( new hdf5IOFile_Type() );
    else
#endif
        M_exporter.reset( new ensightIOFile_Type() );

    M_exporter->setDataFromGetPot( dataFile );
    M_exporter->setPrefix( "Step_" + number2string( multiscaleProblemStep ) + "_Model_" + number2string( M_ID ) );
    M_exporter->setDirectory( multiscaleProblemFolder );

    //Importer
    const std::string importerType = dataFile( "importer/type", "ensight" );

#ifdef HAVE_HDF5
    if ( !importerType.compare( "hdf5" ) )
        M_importer.reset( new hdf5IOFile_Type() );
    else
#endif
        M_importer.reset( new ensightIOFile_Type() );

    M_importer->setDataFromGetPot( dataFile );
    M_importer->setPrefix( "Step_" + number2string( multiscaleProblemStep - 1 ) + "_Model_" + number2string( M_ID ) );
    M_importer->setDirectory( multiscaleProblemFolder );
}

void
MultiscaleModelFluid3D::setupMesh()
{
    //Read fluid mesh from file
    boost::shared_ptr< mesh_Type > fluidMesh( new mesh_Type );
    readMesh( *fluidMesh, *M_dataMesh );

    //Transform mesh
    fluidMesh->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    //Partition mesh
    M_mesh.reset( new partitionMesh_Type( fluidMesh, M_comm ) );
}

void
MultiscaleModelFluid3D::setupFEspace()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::SetupFEspace() \n";
#endif

    //Velocity FE Space
    const RefFE* u_refFE;
    const QuadRule* u_qR;
    const QuadRule* u_bdQr;

    if ( M_data->uOrder().compare( "P2" ) == 0 )
    {
        u_refFE = &feTetraP2;
        u_qR = &quadRuleTetra15pt; // DoE 5
        u_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else if ( M_data->uOrder().compare( "P1" ) == 0 )
    {
        u_refFE = &feTetraP1;
        u_qR = &quadRuleTetra4pt; // DoE 2
        u_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else if ( M_data->uOrder().compare( "P1Bubble" ) == 0 )
    {
        u_refFE = &feTetraP1bubble;
        u_qR = &quadRuleTetra64pt; // DoE 2
        u_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else
    {
        if ( M_displayer->isLeader() )
            std::cout << M_data->uOrder() << " Velocity FE not implemented yet." << std::endl;
        exit( EXIT_FAILURE );
    }

    //Pressure FE Space
    const RefFE* p_refFE;
    const QuadRule* p_qR;
    const QuadRule* p_bdQr;

    if ( M_data->pOrder().compare( "P2" ) == 0 )
    {
        p_refFE = &feTetraP2;
        p_qR = u_qR;
        p_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else if ( M_data->pOrder().compare( "P1" ) == 0 )
    {
        p_refFE = &feTetraP1;
        p_qR = u_qR;
        p_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else
    {
        if ( M_displayer->isLeader() )
            std::cout << M_data->pOrder() << " pressure FE not implemented yet." << std::endl;
        exit( EXIT_FAILURE );
    }

    M_uFESpace.reset( new FESpace_Type( *M_mesh, *u_refFE, *u_qR, *u_bdQr, 3, M_comm ) );
    M_pFESpace.reset( new FESpace_Type( *M_mesh, *p_refFE, *p_qR, *p_bdQr, 1, M_comm ) );
}

void
MultiscaleModelFluid3D::setupDOF()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::SetupDOF \n";
#endif

    M_lmDOF = M_bc->handler()->getNumberBCWithType( Flux );

    //M_uDOF = M_uFESpace->map().getMap(Unique)->NumGlobalElements();
    //M_pFESpace->dof().numTotalDof() = M_pFESpace->map().getMap(Unique)->NumGlobalElements();
}

void
MultiscaleModelFluid3D::setupBCOffset( const bcPtr_Type& BC )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::SetupBCOffset( BC ) \n";
#endif

    UInt offset = M_uFESpace->map().getMap( Unique )->NumGlobalElements() + M_pFESpace->map().getMap( Unique )->NumGlobalElements();

    std::vector< BCName > FluxVector = BC->getBCWithType( Flux );
    for ( UInt i = 0; i < M_lmDOF; ++i )
        BC->setOffset( FluxVector[i], offset + i );
}

void
MultiscaleModelFluid3D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::InitializeSolution() \n";
#endif

    if ( multiscaleProblemStep > 0 )
    {
        M_importer->setMeshProcId( M_mesh->mesh(), M_comm->MyPID() );

        M_importer->addVariable( ExporterData::Vector, "Fluid Velocity", M_solution, static_cast <UInt> ( 0 ),            M_uFESpace->dof().numTotalDof() );
        M_importer->addVariable( ExporterData::Scalar, "Fluid Pressure", M_solution, 3 * M_uFESpace->dof().numTotalDof(), M_pFESpace->dof().numTotalDof());

        // Import
        M_exporter->setStartIndex( M_importer->importFromTime( M_data->dataTime()->getInitialTime() ) + 1 );
    }
    else
        *M_solution = 0.0;

    M_fluid->initialize( *M_solution );
}

void
MultiscaleModelFluid3D::imposePerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::ImposePerturbation() \n";
#endif

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            M_linearBC->GetBCWithFlag( ( *i )->flag( ( *i )->modelGlobalToLocalID( M_ID ) ) ).setBCFunction( M_bcBaseDeltaOne );

            break;
        }
}

void
MultiscaleModelFluid3D::resetPerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8120 ) << "MultiscaleModelFluid3D::ResetPerturbation() \n";
#endif

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            M_linearBC->GetBCWithFlag( ( *i )->flag( ( *i )->modelGlobalToLocalID( M_ID ) ) ).setBCFunction( M_bcBaseDeltaZero );

            break;
        }
}

} // Namespace multiscale
} // Namespace LifeV
