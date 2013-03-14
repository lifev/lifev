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
 *  @brief File containing the Multiscale Model Fluid3D
 *
 *  @date 12-03-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/solver/MultiscaleModelFluid3D.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModelFluid3D::MultiscaleModelFluid3D() :
    multiscaleModel_Type           (),
    MultiscaleInterface            (),
    M_exporter                     (),
    M_importer                     (),
    M_fileName                     (),
    M_fluid                        (),
    M_bc                           ( new bcInterface_Type() ),
    M_bdf                          (),
    M_data                         ( new data_Type() ),
    M_meshData                     ( new MeshData() ),
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
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::MultiscaleModelFluid3D() \n";
#endif

    M_type = Fluid3D;
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModelFluid3D::setupData ( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::setupData( fileName ) \n";
#endif

    multiscaleModel_Type::setupData ( fileName );
    M_fileName = fileName;

    GetPot dataFile ( fileName );

    // Load data
    M_data->setup ( dataFile );
    if ( M_globalData.get() )
    {
        setupGlobalData ( fileName );
    }

    // Mesh data
    M_meshData->setup ( dataFile, "fluid/space_discretization" );

    // Parameters for the NS Iterations
    M_subiterationsMaximumNumber = dataFile ( "fluid/miscellaneous/SubITMax", 0 );
    M_tolerance                  = dataFile ( "fluid/miscellaneous/Tolerance", 1.e-6 );

    M_generalizedAitken.setDefaultOmega (     dataFile ( "fluid/miscellaneous/Omega",        1.e-3 ) );
    M_generalizedAitken.setOmegaMin (         dataFile ( "fluid/miscellaneous/range",        M_generalizedAitken.defaultOmegaFluid() / 1024, 0 ) );
    M_generalizedAitken.setOmegaMax (         dataFile ( "fluid/miscellaneous/range",        M_generalizedAitken.defaultOmegaFluid() * 1024, 1 ) );
    M_generalizedAitken.useDefaultOmega (     dataFile ( "fluid/miscellaneous/fixedOmega",   false ) );
    M_generalizedAitken.setMinimizationType ( dataFile ( "fluid/miscellaneous/inverseOmega", true ) );

    //Boundary Conditions for the problem
    M_bc->fillHandler ( fileName, "fluid" );

    //Setup Exporter & Importer
    setupExporterImporter ( fileName );
}

void
MultiscaleModelFluid3D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::setupModel() \n";
#endif

    //Mesh
    setupMesh();

    //FEspace
    setupFEspace();

    //Add flow rate offset to BC
    M_lmDOF = M_bc->handler()->numberOfBCWithType ( Flux );
    setupBCOffset ( M_bc->handler() );

    //Fluid
    M_fluid.reset ( new fluid_Type ( M_data, *M_uFESpace, *M_pFESpace, M_comm, M_lmDOF ) );
    M_bc->setPhysicalSolver ( M_fluid );

    GetPot dataFile ( M_fileName );
    M_fluid->setUp ( dataFile ); //Remove Preconditioner and Solver if possible!

    //Fluid MAP
    M_map.reset ( new MapEpetra ( M_fluid->getMap() ) );

    //BDF
    M_bdf.reset ( new bdf_Type);
    M_bdf->setup (M_data->dataTimeAdvance()->orderBDF() );

    //Problem coefficients
    M_beta.reset ( new fluidVector_Type ( M_map ) );
    M_rhs.reset ( new fluidVector_Type ( M_map ) );

    //Post-processing
    M_exporter->setMeshProcId ( M_mesh->meshPartition(), M_comm->MyPID() );

    M_solution.reset ( new fluidVector_Type ( *M_fluid->solution(), M_exporter->mapType() ) );
    if ( M_exporter->mapType() == Unique )
    {
        M_solution->setCombineMode ( Zero );
    }

    M_exporter->addVariable ( IOData_Type::VectorField, "Velocity (fluid)", M_uFESpace, M_solution, static_cast<UInt> ( 0 ) );
    M_exporter->addVariable ( IOData_Type::ScalarField, "Pressure (fluid)", M_pFESpace, M_solution, static_cast<UInt> ( 3 * M_uFESpace->dof().numTotalDof() ) );

    //Setup linear model
    setupLinearModel();

    //Setup solution
    initializeSolution();
}

void
MultiscaleModelFluid3D::buildModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::buildModel() \n";
#endif

    // Display data
    //    if ( M_comm->MyPID() == 0 )
    //        M_data->showMe();

    //Build constant matrices
    M_fluid->buildSystem();

    //Initialize BDF
    M_bdf->bdfVelocity().setInitialCondition ( *M_fluid->solution() );

    //Define problem coefficients
    if ( M_data->isStokes() )
    {
        M_alpha  = 0.0;
        *M_beta  = *M_fluid->solution(); //It is a stationary Navier-Stokes
        *M_rhs  *= 0.0;
    }
    else
    {
        M_alpha = M_bdf->bdfVelocity().coefficientFirstDerivative ( 0 ) / M_data->dataTime()->timeStep();
        M_bdf->bdfVelocity().extrapolation (*M_beta);

        M_bdf->bdfVelocity().updateRHSContribution ( M_data->dataTime()->timeStep() );
        *M_rhs  = M_fluid->matrixMass() * M_bdf->bdfVelocity().rhsContributionFirstDerivative();
    }

    //Set problem coefficients
    M_fluid->updateSystem ( M_alpha, *M_beta, *M_rhs );

    //Update operator BC
    M_bc->updatePhysicalSolverVariables();
}

void
MultiscaleModelFluid3D::updateModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::updateModel() \n";
#endif

    //Update BDF
    M_bdf->bdfVelocity().shiftRight ( *M_fluid->solution() );

    //Update problem coefficients
    M_alpha = M_bdf->bdfVelocity().coefficientFirstDerivative ( 0 ) / M_data->dataTime()->timeStep();
    M_bdf->bdfVelocity().extrapolation (*M_beta);

    M_bdf->bdfVelocity().updateRHSContribution ( M_data->dataTime()->timeStep() );
    *M_rhs  = M_fluid->matrixMass() * M_bdf->bdfVelocity().rhsContributionFirstDerivative();

    //Set problem coefficients
    M_fluid->updateSystem ( M_alpha, *M_beta, *M_rhs );

    //Update operator BC
    M_bc->updatePhysicalSolverVariables();

    //Recompute preconditioner
    M_fluid->resetPreconditioner ( true );

    //Linear system need to be updated
    M_updateLinearModel = true;
}

void
MultiscaleModelFluid3D::solveModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::solveModel() \n";
#endif

    // Solve the problem
    displayModelStatus ( "Solve" );
    M_fluid->iterate ( *M_bc->handler() );

    // Non linear convective term with Aitken
    if ( M_subiterationsMaximumNumber > 0 )
    {
        Real residual = ( *M_beta - *M_fluid->solution() ).norm2(); // residual is computed on the whole solution vector;

        if ( M_comm->MyPID() == 0 )
        {
            std::cout << "  F-  Residual:                                " << residual << std::endl;
        }

        M_generalizedAitken.restart();
        for ( UInt subIT = 1; subIT <= M_subiterationsMaximumNumber; ++subIT )
        {
            // Verify tolerance
            if ( residual <= M_tolerance )
            {
                break;
            }

            *M_beta += M_generalizedAitken.computeDeltaLambdaScalar ( *M_beta, *M_beta - *M_fluid->solution() );

            //Linear model need to be updated!
            M_fluid->updateSystem ( M_alpha, *M_beta, *M_rhs );
            M_bc->updatePhysicalSolverVariables();
            M_updateLinearModel = true;

            //Solve system
            M_fluid->iterate ( *M_bc->handler() );

            // Check the new residual
            residual = ( *M_beta - *M_fluid->solution() ).norm2(); // residual is computed on the whole solution vector

            // Display subiteration information
            if ( M_comm->MyPID() == 0 )
            {
                std::cout << "  F-  Sub-iteration n.:                        " << subIT << std::endl;
                std::cout << "  F-  Residual:                                " << residual << std::endl;
            }
        }
    }
}

void
MultiscaleModelFluid3D::updateSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::updateSolution() \n";
#endif

    *M_solution = *M_fluid->solution();
}

void
MultiscaleModelFluid3D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::saveSolution() \n";
#endif

    //Post-processing
    M_exporter->postProcess ( M_data->dataTime()->time() );

#ifdef HAVE_HDF5
    if ( M_data->dataTime()->isLastTimeStep() )
    {
        M_exporter->closeFile();
    }
#endif

}

void
MultiscaleModelFluid3D::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        multiscaleModel_Type::showMe();

        std::cout << "Velocity FE order   = " << M_data->uOrder() << std::endl
                  << "Pressure FE order   = " << M_data->pOrder() << std::endl << std::endl;

        std::cout << "Velocity DOF        = " << 3 * M_uFESpace->dof().numTotalDof() << std::endl
                  << "Pressure DOF        = " << M_pFESpace->dof().numTotalDof() << std::endl
                  << "lmDOF               = " << M_lmDOF << std::endl << std::endl;

        std::cout << "Fluid mesh maxH     = " << MeshUtility::MeshStatistics::computeSize ( *M_mesh->meshPartition() ).maxH << std::endl
                  << "Fluid mesh meanH    = " << MeshUtility::MeshStatistics::computeSize ( *M_mesh->meshPartition() ).meanH << std::endl << std::endl;

        std::cout << "NS SubITMax         = " << M_subiterationsMaximumNumber << std::endl
                  << "NS Tolerance        = " << M_tolerance << std::endl << std::endl << std::endl << std::endl;
    }
}

Real
MultiscaleModelFluid3D::checkSolution() const
{
    return M_solution->norm2();
}

// ===================================================
// MultiscaleInterface Methods
// ===================================================
void
MultiscaleModelFluid3D::imposeBoundaryFlowRate ( const multiscaleID_Type& boundaryID, const function_Type& function )
{
    BCFunctionBase base;
    base.setFunction ( function );

    M_bc->handler()->addBC ( "CouplingFlowRate_Model_" + number2string ( M_ID ) + "_BoundaryID_" + number2string ( boundaryID ), boundaryFlag ( boundaryID ), Flux, Full, base, 3 );
}

void
MultiscaleModelFluid3D::imposeBoundaryMeanNormalStress ( const multiscaleID_Type& boundaryID, const function_Type& function )
{
    BCFunctionBase base;
    base.setFunction ( function );

    M_bc->handler()->addBC ( "CouplingStress_Model_" + number2string ( M_ID ) + "_BoundaryID_" + number2string ( boundaryID ), boundaryFlag ( boundaryID ), Natural, Normal, base );
}

Real
MultiscaleModelFluid3D::boundaryDeltaFlowRate ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem )
{
    solveLinearModel ( solveLinearSystem );

    return M_fluid->linearFlux ( boundaryFlag ( boundaryID ) );
}

Real
MultiscaleModelFluid3D::boundaryDeltaMeanNormalStress ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem )
{
    solveLinearModel ( solveLinearSystem );

    return M_fluid->linearMeanNormalStress ( boundaryFlag ( boundaryID ), *M_linearBC );
}

Real
MultiscaleModelFluid3D::boundaryDeltaMeanTotalNormalStress ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem )
{
    solveLinearModel ( solveLinearSystem );

    return M_fluid->linearMeanTotalNormalStress ( boundaryFlag ( boundaryID ), *M_linearBC );
}

// ===================================================
// Set Methods
// ===================================================
void
MultiscaleModelFluid3D::setSolution ( const fluidVectorPtr_Type& solution )
{
    M_solution = solution;

    M_fluid->initialize ( *M_solution );
}

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleModelFluid3D::setupGlobalData ( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::setupGlobalData( fileName ) \n";
#endif

    GetPot dataFile ( fileName );

    //Global data time
    M_data->setTimeData ( M_globalData->dataTime() );

    //Global physical quantities
    if ( !dataFile.checkVariable ( "fluid/physics/density" ) )
    {
        M_data->setDensity ( M_globalData->fluidDensity() );
    }
    if ( !dataFile.checkVariable ( "fluid/physics/viscosity" ) )
    {
        M_data->setViscosity ( M_globalData->fluidViscosity() );
    }
}

void
MultiscaleModelFluid3D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::initializeSolution() \n";
#endif

    if ( multiscaleProblemStep > 0 )
    {
        M_importer->setMeshProcId ( M_mesh->meshPartition(), M_comm->MyPID() );

        M_importer->addVariable ( IOData_Type::VectorField, "Velocity (fluid)", M_uFESpace, M_solution, static_cast <UInt> ( 0 ) );
        M_importer->addVariable ( IOData_Type::ScalarField, "Pressure (fluid)", M_pFESpace, M_solution, static_cast <UInt> ( 3 * M_uFESpace->dof().numTotalDof() ) );

        // Import
        M_exporter->setTimeIndex ( M_importer->importFromTime ( M_data->dataTime()->initialTime() ) + 1 );

#ifdef HAVE_HDF5
        M_importer->closeFile();
#endif
    }
    else
    {
        *M_solution = 0.0;
    }

    M_fluid->initialize ( *M_solution );
}

void
MultiscaleModelFluid3D::setupExporterImporter ( const std::string& fileName )
{
    GetPot dataFile ( fileName );

    //Exporter
    const std::string exporterType = dataFile ( "exporter/type", "ensight" );

#ifdef HAVE_HDF5
    if ( !exporterType.compare ( "hdf5" ) )
    {
        M_exporter.reset ( new hdf5IOFile_Type() );
    }
    else
#endif
        M_exporter.reset ( new ensightIOFile_Type() );

    M_exporter->setDataFromGetPot ( dataFile );
    M_exporter->setPrefix ( multiscaleProblemPrefix + "_Model_" + number2string ( M_ID ) + "_" + number2string ( multiscaleProblemStep ) );
    M_exporter->setPostDir ( multiscaleProblemFolder );

    //Importer
    const std::string importerType = dataFile ( "importer/type", "ensight" );

#ifdef HAVE_HDF5
    if ( !importerType.compare ( "hdf5" ) )
    {
        M_importer.reset ( new hdf5IOFile_Type() );
    }
    else
#endif
        M_importer.reset ( new ensightIOFile_Type() );

    M_importer->setDataFromGetPot ( dataFile );
    M_importer->setPrefix ( multiscaleProblemPrefix + "_Model_" + number2string ( M_ID ) + "_" + number2string ( multiscaleProblemStep - 1 ) );
    M_importer->setPostDir ( multiscaleProblemFolder );
}

void
MultiscaleModelFluid3D::setupMesh()
{
    //Read fluid mesh from file
    boost::shared_ptr< mesh_Type > fluidMesh ( new mesh_Type ( M_comm ) );
    readMesh ( *fluidMesh, *M_meshData );

    //Transform mesh
    fluidMesh->meshTransformer().transformMesh ( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    //Partition mesh
    M_mesh.reset ( new MeshPartitioner_Type ( fluidMesh, M_comm ) );
}

void
MultiscaleModelFluid3D::setupFEspace()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::setupFEspace() \n";
#endif

    //Velocity FE Space
    const ReferenceFE* u_refFE;
    const QuadratureRule* u_qR;
    const QuadratureRule* u_bdQr;

    if ( M_data->uOrder().compare ( "P2" ) == 0 )
    {
        u_refFE = &feTetraP2;
        u_qR = &quadRuleTetra15pt; // DoE 5
        u_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else if ( M_data->uOrder().compare ( "P1" ) == 0 )
    {
        u_refFE = &feTetraP1;
        u_qR = &quadRuleTetra4pt;  // DoE 2
        u_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else if ( M_data->uOrder().compare ( "P1Bubble" ) == 0 )
    {
        u_refFE = &feTetraP1bubble;
        u_qR = &quadRuleTetra64pt; // DoE 2
        u_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else
    {
        if ( M_comm->MyPID() == 0 )
        {
            std::cout << M_data->uOrder() << " Velocity FE not implemented yet." << std::endl;
        }
        exit ( EXIT_FAILURE );
    }

    //Pressure FE Space
    const ReferenceFE* p_refFE;
    const QuadratureRule* p_qR;
    const QuadratureRule* p_bdQr;

    if ( M_data->pOrder().compare ( "P2" ) == 0 )
    {
        p_refFE = &feTetraP2;
        p_qR = u_qR;
        p_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else if ( M_data->pOrder().compare ( "P1" ) == 0 )
    {
        p_refFE = &feTetraP1;
        p_qR = u_qR;
        p_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else
    {
        if ( M_comm->MyPID() == 0 )
        {
            std::cout << M_data->pOrder() << " pressure FE not implemented yet." << std::endl;
        }
        exit ( EXIT_FAILURE );
    }

    M_uFESpace.reset ( new FESpace_Type ( M_mesh->meshPartition(), *u_refFE, *u_qR, *u_bdQr, 3, M_comm ) );
    M_pFESpace.reset ( new FESpace_Type ( M_mesh->meshPartition(), *p_refFE, *p_qR, *p_bdQr, 1, M_comm ) );
}

void
MultiscaleModelFluid3D::setupDOF()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::setupDOF \n";
#endif

    M_lmDOF = M_bc->handler()->numberOfBCWithType ( Flux );
}

void
MultiscaleModelFluid3D::setupBCOffset ( const bcPtr_Type& bc )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::setupBCOffset( BC ) \n";
#endif

    UInt offset = M_uFESpace->map().map ( Unique )->NumGlobalElements() + M_pFESpace->map().map ( Unique )->NumGlobalElements();

    std::vector< bcName_Type > fluxVector = bc->findAllBCWithType ( Flux );
    for ( UInt i = 0; i < M_lmDOF; ++i )
    {
        bc->setOffset ( fluxVector[i], offset + i );
    }
}

void
MultiscaleModelFluid3D::setupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::setupLinearModel( ) \n";
#endif

    // The linear BCHandler is a copy of the original BCHandler with all BCFunctions giving zero
    M_linearBC.reset ( new bc_Type ( *M_bc->handler() ) );

    // Set all the BCFunctions to zero
    BCFunctionBase bcBaseDeltaZero;
    bcBaseDeltaZero.setFunction ( boost::bind ( &MultiscaleModelFluid3D::bcFunctionDeltaZero, this, _1, _2, _3, _4, _5 ) );

    for ( bc_Type::bcBaseIterator_Type i = M_linearBC->begin() ; i != M_linearBC->end() ; ++i )
    {
        i->setBCFunction ( bcBaseDeltaZero );
    }
}

void
MultiscaleModelFluid3D::updateLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::updateLinearModel() \n";
#endif

    //Create an empty vector
    fluidVector_Type vectorZero ( *M_solution );
    vectorZero = 0.0;

    //updateLinearModel TODO REMOVE ?
    M_fluid->updateLinearSystem ( M_fluid->matrixNoBC(), M_alpha, *M_beta, *M_fluid->solution(),
                                  vectorZero, vectorZero, vectorZero, vectorZero );

    //Linear System Updated
    M_updateLinearModel = false;
}

void
MultiscaleModelFluid3D::solveLinearModel ( bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::solveLinearModel() \n";
#endif

    if ( !solveLinearSystem )
    {
        return;
    }

    imposePerturbation();

    if ( M_updateLinearModel )
    {
        updateLinearModel();
    }

    //Solve the linear problem
    displayModelStatus ( "Solve linear" );
    M_fluid->solveLinearSystem ( *M_linearBC );

    resetPerturbation();

    //This flag avoid recomputation of the same system
    solveLinearSystem = false;
}

void
MultiscaleModelFluid3D::imposePerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::imposePerturbation() \n";
#endif

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            BCFunctionBase bcBaseDeltaOne;
            bcBaseDeltaOne.setFunction ( boost::bind ( &MultiscaleModelFluid3D::bcFunctionDeltaOne, this, _1, _2, _3, _4, _5 ) );

            M_linearBC->findBCWithFlag ( boundaryFlag ( ( *i )->boundaryID ( ( *i )->modelGlobalToLocalID ( M_ID ) ) ) ).setBCFunction ( bcBaseDeltaOne );

            break;
        }
}

void
MultiscaleModelFluid3D::resetPerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8120 ) << "MultiscaleModelFluid3D::resetPerturbation() \n";
#endif

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            BCFunctionBase bcBaseDeltaZero;
            bcBaseDeltaZero.setFunction ( boost::bind ( &MultiscaleModelFluid3D::bcFunctionDeltaZero, this, _1, _2, _3, _4, _5 ) );

            M_linearBC->findBCWithFlag ( boundaryFlag ( ( *i )->boundaryID ( ( *i )->modelGlobalToLocalID ( M_ID ) ) ) ).setBCFunction ( bcBaseDeltaZero );

            break;
        }
}

} // Namespace multiscale
} // Namespace LifeV
