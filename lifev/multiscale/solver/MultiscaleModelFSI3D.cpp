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

#include <lifev/multiscale/solver/MultiscaleModelFSI3D.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModelFSI3D::MultiscaleModelFSI3D() :
    multiscaleModel_Type           (),
    MultiscaleInterface            (),
    M_FSIoperator                  (),
    M_data                         ( new data_Type() ),
    M_exporterFluid                (),
    M_exporterSolid                (),
    M_importerFluid                (),
    M_importerSolid                (),
#ifndef FSI_WITH_EXTERNALPRESSURE
    M_boundaryStressFunctions      (),
    M_externalPressureScalar       (),
#endif
#ifndef FSI_WITHOUT_VELOCITYPROFILE
    M_boundaryFlowRateFunctions    (),
    M_boundaryFlowRateType         (),
#endif
#ifdef FSI_WITH_BOUNDARYAREA
    M_boundaryAreaFunctions        (),
    M_boundaryFlagsArea            (),
#endif
    M_fluidVelocity                (),
    M_fluidPressure                (),
    M_fluidDisplacement            (),
    M_solidVelocity                (),
    M_solidDisplacement            (),
    M_nonLinearRichardsonIteration ( 0 ),
    M_fluidBC                      ( new bcInterface_Type() ),
    M_solidBC                      ( new bcInterface_Type() ),
    M_harmonicExtensionBC          ( new bcInterface_Type() ),
    M_linearFluidBC                (),
    M_linearSolidBC                (),
    M_linearRHS                    (),
    M_linearSolution               ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::MultiscaleModelFSI3D() \n";
#endif

    M_type = FSI3D;

    FSIOperator_Type::factory_Type::instance().registerProduct ( "monolithicGE", &createFSIMonolithicGE );
    FSIOperator_Type::factory_Type::instance().registerProduct ( "monolithicGI", &createFSIMonolithicGI );

    BlockPrecFactory::instance().registerProduct ("AdditiveSchwarz",   &MonolithicBlockMatrix::createAdditiveSchwarz) ;
    BlockPrecFactory::instance().registerProduct ("AdditiveSchwarzRN", &MonolithicBlockMatrixRN::createAdditiveSchwarzRN ) ;
    BlockPrecFactory::instance().registerProduct ("ComposedDN",        &MonolithicBlockComposedDN::createComposedDN ) ;
    BlockPrecFactory::instance().registerProduct ("ComposedDN2",       &MonolithicBlockComposedDN::createComposedDN2 );
    BlockPrecFactory::instance().registerProduct ("ComposedDNND",      &MonolithicBlockComposedDNND::createComposedDNND);

    MonolithicBlockMatrix::Factory_Type::instance().registerProduct ("AdditiveSchwarz",   &MonolithicBlockMatrix::createAdditiveSchwarz ) ;
    MonolithicBlockMatrix::Factory_Type::instance().registerProduct ("AdditiveSchwarzRN", &MonolithicBlockMatrixRN::createAdditiveSchwarzRN ) ;

    FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct ( "linearVenantKirchhoff",    &FSIOperator::createVenantKirchhoffLinear );
    FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct ( "nonLinearVenantKirchhoff", &FSIOperator::createVenantKirchhoffNonLinear );
    FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct ( "exponential",              &FSIOperator::createExponentialMaterialNonLinear );
    FSIOperator::solid_Type::material_Type::StructureMaterialFactory::instance().registerProduct ( "neoHookean",               &FSIOperator::createNeoHookeanMaterialNonLinear );
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModelFSI3D::setupData ( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::setupData( fileName ) \n";
#endif

    multiscaleModel_Type::setupData ( fileName );

    GetPot dataFile ( fileName );

    // Load data
    M_data->setup ( dataFile );
    if ( M_globalData.get() )
    {
        setupGlobalData ( fileName );
    }

    // Create FSI
    M_FSIoperator = FSIOperatorPtr_Type ( FSIOperator_Type::factory_Type::instance().createObject ( M_data->method() ) );

    // Setup Communicator
    setupCommunicator();

    // Set data
    M_FSIoperator->setData ( M_data );
    M_FSIoperator->setDataFile ( dataFile );

    // Setup Boundary Conditions from file
    setupBC ( fileName );

    // Setup exporter
    if ( M_FSIoperator->isFluid() )
    {
        setupExporter ( M_exporterFluid, dataFile );
        setupImporter ( M_importerFluid, dataFile );
    }
    if ( M_FSIoperator->isSolid() )
    {
        setupExporter ( M_exporterSolid, dataFile, "_Solid" );
        setupImporter ( M_importerSolid, dataFile, "_Solid" );
    }

#ifndef FSI_WITHOUT_VELOCITYPROFILE
    std::map< std::string, FSI3DBoundaryFlowRate_Type > boundaryFlowRateMap;
    boundaryFlowRateMap["Weak"]     = Weak;
    boundaryFlowRateMap["Semiweak"] = Semiweak;
    boundaryFlowRateMap["Strong"]   = Strong;

    M_boundaryFlowRateType.reserve ( M_boundaryFlags.size() );
    for ( UInt j ( 0 ); j < M_boundaryFlags.size(); ++j )
    {
        M_boundaryFlowRateType[j] = boundaryFlowRateMap[dataFile ( "Multiscale/flowRateCouplingType", "Weak", j )];
    }
#endif

#ifdef FSI_WITH_BOUNDARYAREA
    M_boundaryFlagsArea.reserve ( M_boundaryFlags.size() );
    M_boundaryFlagsAreaPerturbed.reserve ( M_boundaryFlags.size() );
    for ( UInt j ( 0 ); j < M_boundaryFlags.size(); ++j )
    {
        M_boundaryFlagsArea.push_back ( dataFile ( "Multiscale/couplingAreaFlags", 0, j ) );
        M_boundaryFlagsAreaPerturbed.push_back ( false );
    }
#endif
}

void
MultiscaleModelFSI3D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::setupProblem() \n";
#endif

    // Mesh transformation (before partitioning, ideally should be done after for scalability)
    //M_FSIoperator->fluidMesh().transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );
    M_FSIoperator->solidMesh().meshTransformer().transformMesh ( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    // Mesh partitioning
    M_FSIoperator->partitionMeshes();

    // Mesh transformation (after partitioning - not working for solid)
    M_FSIoperator->fluidLocalMesh().meshTransformer().transformMesh ( M_geometryScale, M_geometryRotate, M_geometryTranslate );
    //M_FSIoperator->solidLocalMesh().meshTransformer()->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    // Setup FEspace & DOF
    M_FSIoperator->setupFEspace();
    M_FSIoperator->setupDOF();

    // Setup FSI Interface Boundary Conditions (by giving the operator to BCInterface)
    M_fluidBC->setPhysicalSolver ( M_FSIoperator );
    M_solidBC->setPhysicalSolver ( M_FSIoperator );
    M_harmonicExtensionBC->setPhysicalSolver ( M_FSIoperator );

    // Setup Fluid & Solid solver
    M_FSIoperator->setupFluidSolid ( M_FSIoperator->imposedFluxes() );
    M_FSIoperator->setupSystem();

#ifdef FSI_WITH_BOUNDARYAREA
    // Setup the boundary area functions
    for ( boundaryAreaFunctionsContainerIterator_Type i = M_boundaryAreaFunctions.begin(); i < M_boundaryAreaFunctions.end(); ++i )
    {
        ( *i )->setup();
    }
#endif

    // Setup Exporters
    if ( M_FSIoperator->isFluid() )
    {
        setExporterFluid ( M_exporterFluid );
    }

    if ( M_FSIoperator->isSolid() )
    {
        setExporterSolid ( M_exporterSolid );
    }

    //Setup solution
    initializeSolution();

    //Setup linear model
    setupLinearModel();
}

void
MultiscaleModelFSI3D::buildModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::buildModel() \n";
#endif

    // Display data
    //    if ( M_comm->MyPID() == 0 )
    //        M_data->showMe();

    M_FSIoperator->buildSystem();
    M_FSIoperator->updateSystem();

    // Update BCInterface solver variables
    updateBC();
}

void
MultiscaleModelFSI3D::updateModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::updateModel() \n";
#endif

    // Update System
    M_FSIoperator->updateSystem();

    // Update BCInterface solver variables
    updateBC();

    // TODO This is a workaround of Paolo Crosetto to make it possible to perform subiterations
    // in the multiscale when using 3D FSI models. In the future this should be replaced with
    // a more proper implementation.
    M_FSIoperator->precPtrView()->setRecompute ( 1, true );
    M_nonLinearRichardsonIteration = 0;
}

void
MultiscaleModelFSI3D::solveModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::solveModel() \n";
#endif

    displayModelStatus ( "Solve" );

    // TODO This is a workaround of Paolo Crosetto to make it possible to perform subiterations
    // in the multiscale when using 3D FSI models. In the future this should be replaced with
    // a more proper implementation.
    if ( M_nonLinearRichardsonIteration != 0 )
    {
        M_FSIoperator->resetRHS();
        M_FSIoperator->updateRHS();
        M_FSIoperator->applyBoundaryConditions();
    }

    // Non-linear Richardson solver
    UInt maxSubIterationNumber = M_data->maxSubIterationNumber();
    M_FSIoperator->extrapolation ( *M_stateVariable );

    NonLinearRichardson ( *M_stateVariable, *M_FSIoperator,
                          M_data->absoluteTolerance(), M_data->relativeTolerance(),
                          maxSubIterationNumber, M_data->errorTolerance(),
                          M_data->NonLinearLineSearch(),
                          M_nonLinearRichardsonIteration,
                          1
                        );

    // TODO This is a workaround of Paolo Crosetto to make it possible to perform subiterations
    // in the multiscale when using 3D FSI models. In the future this should be replaced with
    // a more proper implementation.
    M_FSIoperator->precPtrView()->setRecompute ( 1, false );
    M_nonLinearRichardsonIteration = 1;

#ifdef FSI_WITH_BOUNDARYAREA
    for ( UInt j ( 0 ); j < M_boundaryFlagsAreaPerturbed.size(); ++j )
    {
        M_boundaryFlagsAreaPerturbed[j] = false;
    }
#endif
}

void
MultiscaleModelFSI3D::updateSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::updateSolution() \n";
#endif

    M_FSIoperator->updateSolution ( *M_stateVariable );
}

void
MultiscaleModelFSI3D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::saveSolution() \n";
#endif

    exportFluidSolution();
    if ( M_FSIoperator->isFluid() )
    {
        M_exporterFluid->postProcess ( M_data->dataFluid()->dataTime()->time() );
    }

    exportSolidSolution();
    if ( M_FSIoperator->isSolid() )
    {
        M_exporterSolid->postProcess ( M_data->dataSolid()->dataTime()->time() );
    }

#ifdef HAVE_HDF5
    if ( M_data->dataFluid()->dataTime()->isLastTimeStep() )
    {
        if ( M_FSIoperator->isFluid() )
        {
            M_exporterFluid->closeFile();
        }
        if ( M_FSIoperator->isSolid() )
        {
            M_exporterSolid->closeFile();
        }
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
                  << "Structure FE order  = " << M_FSIoperator->dataSolid()->order() << std::endl << std::endl;

        std::cout << "Velocity DOF        = " << 3 * M_FSIoperator->uFESpace().dof().numTotalDof() << std::endl
                  << "Pressure DOF        = " << M_FSIoperator->pFESpace().dof().numTotalDof() << std::endl
                  << "Harmonic ext. DOF   = " << M_FSIoperator->mmFESpace().dof().numTotalDof() << std::endl
                  << "Structure DOF       = " << M_FSIoperator->dFESpace().dof().numTotalDof() << std::endl << std::endl;

        std::cout << "Fluid mesh maxH     = " << MeshUtility::MeshStatistics::computeSize ( * ( M_FSIoperator->uFESpace().mesh() ) ).maxH << std::endl
                  << "Fluid mesh meanH    = " << MeshUtility::MeshStatistics::computeSize ( * ( M_FSIoperator->uFESpace().mesh() ) ).meanH << std::endl
                  << "Solid mesh maxH     = " << MeshUtility::MeshStatistics::computeSize ( * ( M_FSIoperator->dFESpace().mesh() ) ).maxH << std::endl
                  << "Solid mesh meanH    = " << MeshUtility::MeshStatistics::computeSize ( * ( M_FSIoperator->dFESpace().mesh() ) ).meanH << std::endl << std::endl;
    }
}

Real
MultiscaleModelFSI3D::checkSolution() const
{
    return M_stateVariable->norm2();
}

// ===================================================
// MultiscaleInterface Methods
// ===================================================
void
MultiscaleModelFSI3D::imposeBoundaryFlowRate ( const multiscaleID_Type& boundaryID, const function_Type& function )
{
    BCFunctionBase base;

#ifndef FSI_WITHOUT_VELOCITYPROFILE
    switch ( M_boundaryFlowRateType[boundaryID] )
    {
        case Strong:
        {
            boundaryFlowRateFunctionPtr_Type boundaryFlowRateFunction ( new boundaryFlowRateFunction_Type() );
            boundaryFlowRateFunction->setModel ( this );
            boundaryFlowRateFunction->setFluidFlag ( boundaryFlag ( boundaryID ) );
            boundaryFlowRateFunction->setFunction ( function);
            boundaryFlowRateFunction->setBoundaryFlowRateType ( Strong );

            M_boundaryFlowRateFunctions.push_back ( boundaryFlowRateFunction );

            base.setFunction ( boost::bind ( &FSI3DBoundaryFlowRateFunction::function, M_boundaryFlowRateFunctions.back(), _1, _2, _3, _4, _5 ) );
            M_fluidBC->handler()->addBC ( "CouplingFlowRate_Model_" + number2string ( M_ID ) + "_BoundaryID_" + number2string ( boundaryID ), boundaryFlag ( boundaryID ), EssentialVertices, Full, base, 3 );

            break;
        }
        case Semiweak:
        {
            // TODO: implementation is still required for the switching of the boundary condition.
            // Possible idea: we can use the same function and change just the type of BC inside the BCHandler from Flux to Essential.
            boundaryFlowRateFunctionPtr_Type boundaryFlowRateFunction ( new boundaryFlowRateFunction_Type() );
            boundaryFlowRateFunction->setModel ( this );
            boundaryFlowRateFunction->setFluidFlag ( boundaryFlag ( boundaryID ) );
            boundaryFlowRateFunction->setFunction ( function);
            boundaryFlowRateFunction->setBoundaryFlowRateType ( Semiweak );

            M_boundaryFlowRateFunctions.push_back ( boundaryFlowRateFunction );

            base.setFunction ( boost::bind ( &FSI3DBoundaryFlowRateFunction::function, M_boundaryFlowRateFunctions.back(), _1, _2, _3, _4, _5 ) );
            M_fluidBC->handler()->addBC ( "CouplingFlowRate_Model_" + number2string ( M_ID ) + "_BoundaryID_" + number2string ( boundaryID ), boundaryFlag ( boundaryID ), Flux, Full, base, 3 );

            break;
        }
        case Weak:
        default:

            base.setFunction ( function );
            M_fluidBC->handler()->addBC ( "CouplingFlowRate_Model_" + number2string ( M_ID ) + "_BoundaryID_" + number2string ( boundaryID ), boundaryFlag ( boundaryID ), Flux, Full, base, 3 );

            break;
    }
#else
    base.setFunction ( function );
    M_fluidBC->handler()->addBC ( "CouplingFlowRate_Model_" + number2string ( M_ID ) + "_BoundaryID_" + number2string ( boundaryID ), boundaryFlag ( boundaryID ), Flux, Full, base, 3 );
#endif
}

void
MultiscaleModelFSI3D::imposeBoundaryMeanNormalStress ( const multiscaleID_Type& boundaryID, const function_Type& function )
{
    BCFunctionBase base;
#ifdef FSI_WITH_EXTERNALPRESSURE
    base.setFunction ( function );
#else
    boundaryStressFunctionPtr_Type boundaryStressFunction ( new boundaryStressFunction_Type() );
    boundaryStressFunction->setDelta ( M_data->dataSolid()->externalPressure() );
    boundaryStressFunction->setFunction ( function);

    M_boundaryStressFunctions.push_back ( boundaryStressFunction );

    base.setFunction ( boost::bind ( &FSI3DBoundaryStressFunction::function, M_boundaryStressFunctions.back(), _1, _2, _3, _4, _5 ) );
#endif
    M_fluidBC->handler()->addBC ( "BoundaryStress_Model_" + number2string ( M_ID ) + "_BoundaryID_" + number2string ( boundaryID ), boundaryFlag ( boundaryID ), Natural, Normal, base );
}

void
MultiscaleModelFSI3D::imposeBoundaryArea ( const multiscaleID_Type& boundaryID, const function_Type& function )
{
#ifdef FSI_WITH_BOUNDARYAREA
    boundaryAreaFunctionPtr_Type boundaryAreaFunction ( new boundaryAreaFunction_Type() );
    boundaryAreaFunction->setModel ( this );
    boundaryAreaFunction->setFluidFlag ( boundaryFlag ( boundaryID ) );
    boundaryAreaFunction->setFunction ( function);

    M_boundaryAreaFunctions.push_back ( boundaryAreaFunction );

    BCFunctionBase base;
    base.setFunction ( boost::bind ( &FSI3DBoundaryAreaFunction::function, M_boundaryAreaFunctions.back(), _1, _2, _3, _4, _5 ) );
    M_solidBC->handler()->addBC ( "BoundaryArea_Model_" + number2string ( M_ID ) + "_BoundaryID_" + number2string ( boundaryID ), M_boundaryFlagsArea[boundaryID], EssentialEdges, Full, base, 3 );
#endif
}

Real
MultiscaleModelFSI3D::boundaryMeanNormalStress ( const multiscaleID_Type& boundaryID ) const
{
#ifdef FSI_WITH_EXTERNALPRESSURE
    return M_FSIoperator->fluid().meanNormalStress ( boundaryFlag ( boundaryID ), *M_fluidBC->handler(), *M_stateVariable );
#else
    return M_FSIoperator->fluid().meanNormalStress ( boundaryFlag ( boundaryID ), *M_fluidBC->handler(), *M_stateVariable ) - M_externalPressureScalar;
#endif
}

Real
MultiscaleModelFSI3D::boundaryMeanTotalNormalStress ( const multiscaleID_Type& boundaryID ) const
{
#ifdef FSI_WITH_EXTERNALPRESSURE
    return M_FSIoperator->fluid().meanTotalNormalStress ( boundaryFlag ( boundaryID ), *M_fluidBC->handler(), *M_stateVariable );
#else
    return M_FSIoperator->fluid().meanTotalNormalStress ( boundaryFlag ( boundaryID ), *M_fluidBC->handler(), *M_stateVariable ) - M_externalPressureScalar;
#endif
}

Real
MultiscaleModelFSI3D::boundaryDeltaFlowRate ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem )
{
    solveLinearModel ( solveLinearSystem );

    return M_FSIoperator->fluid().flux ( boundaryFlag ( boundaryID ), *M_linearSolution );
}

Real
MultiscaleModelFSI3D::boundaryDeltaMeanNormalStress ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem )
{
    solveLinearModel ( solveLinearSystem );

    return M_FSIoperator->fluid().linearMeanNormalStress ( boundaryFlag ( boundaryID ), *M_linearFluidBC, *M_linearSolution );
}

Real
MultiscaleModelFSI3D::boundaryDeltaMeanTotalNormalStress ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem )
{
    solveLinearModel ( solveLinearSystem );

    return M_FSIoperator->fluid().linearMeanTotalNormalStress ( boundaryFlag ( boundaryID ), *M_linearFluidBC, *M_stateVariable, *M_linearSolution );
}

// ===================================================
// Get Methods
// ===================================================
Real
MultiscaleModelFSI3D::boundaryPressure ( const multiscaleID_Type& boundaryID ) const
{
#ifdef FSI_WITH_EXTERNALPRESSURE
    return M_FSIoperator->fluid().pressure ( boundaryFlag ( boundaryID ), *M_stateVariable );
#else
    return M_FSIoperator->fluid().pressure ( boundaryFlag ( boundaryID ), *M_stateVariable ) + M_externalPressureScalar;
#endif
}

Real
MultiscaleModelFSI3D::boundaryTotalPressure ( const multiscaleID_Type& boundaryID ) const
{
#ifdef FSI_WITH_EXTERNALPRESSURE
    return M_FSIoperator->fluid().pressure ( boundaryFlag ( boundaryID ), *M_stateVariable )
           + M_FSIoperator->fluid().kineticNormalStress ( boundaryFlag ( boundaryID ), *M_stateVariable );
#else
    return M_FSIoperator->fluid().pressure ( boundaryFlag ( boundaryID ), *M_stateVariable )
           + M_FSIoperator->fluid().kineticNormalStress ( boundaryFlag ( boundaryID ), *M_stateVariable ) + M_externalPressureScalar;
#endif
}

Real
MultiscaleModelFSI3D::externalPressure() const
{
#ifndef FSI_WITH_EXTERNALPRESSURE
    return M_externalPressureScalar;
#else
    return M_data->dataSolid()->externalPressure();
#endif
}

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleModelFSI3D::setupGlobalData ( const std::string& fileName )
{
    GetPot dataFile ( fileName );

    // Global data time
    M_data->dataFluid()->setTimeData ( M_globalData->dataTime() );
    M_data->dataSolid()->setTimeData ( M_globalData->dataTime() );
    M_data->setTimeDataALE ( M_globalData->dataTime() );

    // Fluid global physical quantities
    if ( !dataFile.checkVariable ( "fluid/physics/density" ) )
    {
        M_data->dataFluid()->setDensity ( M_globalData->fluidDensity() );
    }
    if ( !dataFile.checkVariable ( "fluid/physics/viscosity" ) )
    {
        M_data->dataFluid()->setViscosity ( M_globalData->fluidViscosity() );
    }

    // Solid global physical quantities
    if ( !dataFile.checkVariable ( "solid/physics/density" ) )
    {
        M_data->dataSolid()->setDensity ( M_globalData->solidDensity() );
    }
    if ( !dataFile.checkVariable ( "solid/physics/externalPressure" ) )
    {
        M_data->dataSolid()->setExternalPressure ( M_globalData->solidExternalPressure() );
    }

    std::vector< UInt > materialFlags;
    if ( !dataFile.checkVariable ( "solid/physics/material_flag" ) )
    {
        materialFlags.push_back ( 1 );
    }
    else
        for ( UInt i ( 0 ) ; i < dataFile.vector_variable_size ( "solid/physics/material_flag" ) ; ++i )
        {
            materialFlags.push_back ( dataFile ( "solid/physics/material_flag", 1, i ) );
        }

    for ( std::vector< UInt >::const_iterator i = materialFlags.begin(); i != materialFlags.end() ; ++i )
    {
        if ( !dataFile.checkVariable ( "solid/physics/poisson" ) )
        {
            M_data->dataSolid()->setPoisson ( M_globalData->solidPoissonCoefficient(), *i );
        }
        if ( !dataFile.checkVariable ( "solid/physics/young" ) )
        {
            M_data->dataSolid()->setYoung ( M_globalData->solidYoungModulus(), *i );
        }
    }
}

void
MultiscaleModelFSI3D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::initializeSolution() \n";
#endif

#ifndef FSI_WITH_EXTERNALPRESSURE
    // Initialize the external pressure scalar
    M_externalPressureScalar = M_data->dataSolid()->externalPressure();
    M_data->dataSolid()->setExternalPressure ( 0.0 );
#endif

    if ( multiscaleProblemStep > 0 )
    {
        M_importerFluid->setMeshProcId ( M_FSIoperator->uFESpace().mesh(), M_FSIoperator->uFESpace().map().comm().MyPID() );
        M_importerSolid->setMeshProcId ( M_FSIoperator->dFESpace().mesh(), M_FSIoperator->dFESpace().map().comm().MyPID() );

        M_importerFluid->addVariable ( IOData_Type::VectorField, "Velocity (fluid)",     M_FSIoperator->uFESpacePtr(),  M_fluidVelocity,     static_cast <UInt> (0) );
        M_importerFluid->addVariable ( IOData_Type::ScalarField, "Pressure (fluid)",     M_FSIoperator->pFESpacePtr(),  M_fluidPressure,     static_cast <UInt> (0) );
        M_importerFluid->addVariable ( IOData_Type::VectorField, "Displacement (fluid)", M_FSIoperator->mmFESpacePtr(), M_fluidDisplacement, static_cast <UInt> (0) );
        M_importerSolid->addVariable ( IOData_Type::VectorField, "Displacement (solid)", M_FSIoperator->dFESpacePtr(),  M_solidDisplacement, static_cast <UInt> (0) );

        // Read fluid velocity and pressure + solid displacement
        for ( UInt i (0); i < M_FSIoperator->fluidTimeAdvance()->size() ; ++i )
        {
            UInt iterationImported (0);
            iterationImported = M_importerFluid->importFromTime ( M_data->dataFluid()->dataTime()->initialTime() - i * M_data->dataFluid()->dataTime()->timeStep() );
            iterationImported = M_importerSolid->importFromTime ( M_data->dataSolid()->dataTime()->initialTime() - i * M_data->dataSolid()->dataTime()->timeStep() );
            if ( i == 0 )
            {
                M_exporterFluid->setTimeIndex ( iterationImported + 1 );
                M_exporterSolid->setTimeIndex ( iterationImported + 1 );
            }

#ifndef FSI_WITH_EXTERNALPRESSURE
            // Remove external pressure
            *M_fluidPressure -= M_externalPressureScalar;
#endif

            M_FSIoperator->setFluidVectorInStencil ( M_fluidVelocity, M_fluidPressure, i );
            M_FSIoperator->setSolidVectorInStencil ( M_solidDisplacement, i );
        }

        // Read last solid displacement
        M_importerSolid->importFromTime ( M_data->dataSolid()->dataTime()->initialTime() - M_FSIoperator->fluidTimeAdvance()->size() * M_data->dataSolid()->dataTime()->timeStep() );
        M_FSIoperator->setSolidVectorInStencil ( M_solidDisplacement, M_FSIoperator->fluidTimeAdvance()->size() );

        // Read fluid displacement
        for ( UInt i (0); i < M_FSIoperator->ALETimeAdvance()->size() ; ++i )
        {
            M_importerFluid->importFromTime ( M_data->dataFluid()->dataTime()->initialTime() - (i + 1) * M_data->dataFluid()->dataTime()->timeStep() );

            M_FSIoperator->setALEVectorInStencil ( M_fluidDisplacement, i, false );
        }

        // Read first fluid displacement
        M_importerFluid->importFromTime ( M_data->dataSolid()->dataTime()->initialTime() );

        //This is ugly but it's the only way I have figured out at the moment
        if ( M_data->method().compare ("monolithicGI") == 0 )
        {
            //Don't be scared by the ten. The goal of 10 is just to make the first if fail
            M_FSIoperator->setALEVectorInStencil ( M_fluidDisplacement, 10, true );
        }

        //Setting the vector in the stencil
        M_FSIoperator->ALETimeAdvance()->shiftRight ( *M_fluidDisplacement );

        // Finalize restart
        M_FSIoperator->finalizeRestart();

#ifdef HAVE_HDF5
        if ( M_FSIoperator->isFluid() )
        {
            M_importerFluid->closeFile();
        }
        if ( M_FSIoperator->isSolid() )
        {
            M_importerSolid->closeFile();
        }
#endif

    }
    else
    {
        // Initialize containers to zero
        *M_fluidVelocity     *= 0.0;
        *M_fluidPressure      = M_data->dataSolid()->externalPressure();
        *M_solidDisplacement *= 0.0;
        *M_fluidDisplacement *= 0.0;

        for ( UInt i (0); i < M_FSIoperator->fluidTimeAdvance()->size() ; ++i )
        {
            M_FSIoperator->setFluidVectorInStencil ( M_fluidVelocity, M_fluidPressure, i );
            M_FSIoperator->setSolidVectorInStencil ( M_solidDisplacement, i );
        }

        M_FSIoperator->setSolidVectorInStencil ( M_solidDisplacement, M_FSIoperator->fluidTimeAdvance()->size() );

        for ( UInt i (0); i < M_FSIoperator->ALETimeAdvance()->size() ; ++i )
        {
            M_FSIoperator->setALEVectorInStencil ( M_fluidDisplacement, i, false );
        }

        //This is ugly but it's the only way I have figured out at the moment
        if ( M_data->method().compare ("monolithicGI") == 0 )
        {
            //Don't be scared by the ten. The goal of 10 is just to make the first if fail
            M_FSIoperator->setALEVectorInStencil ( M_fluidDisplacement, 10, true );
        }

        //Setting the vector in the stencil
        M_FSIoperator->ALETimeAdvance()->shiftRight ( *M_fluidDisplacement );

        // Finalize restart
        M_FSIoperator->finalizeRestart();
    }

    // Re-initialize fluid and velocity vectors
    exportFluidSolution();
    exportSolidSolution();

    // Initialize state variables
    M_stateVariable.reset ( new vector_Type ( M_FSIoperator->fluidTimeAdvance()->singleElement (0) ) );
}

void
MultiscaleModelFSI3D::setupCommunicator()
{
    M_FSIoperator->setFluid ( true );
    M_FSIoperator->setSolid ( true );

    M_FSIoperator->setFluidLeader ( 0 );
    M_FSIoperator->setSolidLeader ( 0 );

    M_FSIoperator->setComm ( M_comm, M_comm );
}

void
MultiscaleModelFSI3D::setupBC ( const std::string& fileName )
{
    if ( M_FSIoperator->isFluid() )
    {
        M_fluidBC->createHandler();
        M_fluidBC->fillHandler ( fileName, "fluid" );

        M_FSIoperator->setFluidBC ( M_fluidBC->handler() );

        M_harmonicExtensionBC->createHandler();
        M_harmonicExtensionBC->fillHandler ( fileName, "mesh_motion" );

        M_FSIoperator->setHarmonicExtensionBC ( M_harmonicExtensionBC->handler() );
    }

    if ( M_FSIoperator->isSolid() )
    {
        M_solidBC->createHandler();
        M_solidBC->fillHandler ( fileName, "solid" );

        M_FSIoperator->setSolidBC ( M_solidBC->handler() );
    }
}

void
MultiscaleModelFSI3D::updateBC()
{
#ifndef FSI_WITHOUT_VELOCITYPROFILE
    for ( boundaryFlowRateFunctionsContainerIterator_Type i = M_boundaryFlowRateFunctions.begin() ; i != M_boundaryFlowRateFunctions.end() ; ++i )
    {
        ( *i )->updateParameters();
    }
#endif

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
MultiscaleModelFSI3D::setupExporter ( IOFilePtr_Type& exporter, const GetPot& dataFile, const std::string& label )
{
    const std::string exporterType = dataFile ( "exporter/type", "ensight" );
#ifdef HAVE_HDF5
    if ( exporterType.compare ( "hdf5" ) == 0 )
    {
        exporter.reset ( new hdf5IOFile_Type() );
    }
    else
#endif
        exporter.reset ( new ensightIOFile_Type() );

    exporter->setDataFromGetPot ( dataFile );
    exporter->setPrefix ( multiscaleProblemPrefix + "_Model_" + number2string ( M_ID ) +  label + "_" + number2string ( multiscaleProblemStep ) );
    exporter->setPostDir ( multiscaleProblemFolder );
}

void
MultiscaleModelFSI3D::setupImporter ( IOFilePtr_Type& importer, const GetPot& dataFile, const std::string& label )
{
    const std::string importerType = dataFile ( "exporter/type", "ensight" );
#ifdef HAVE_HDF5
    if ( importerType.compare ( "hdf5" ) == 0 )
    {
        importer.reset ( new hdf5IOFile_Type() );
    }
    else
#endif
        importer.reset ( new ensightIOFile_Type() );

    importer->setDataFromGetPot ( dataFile );
    importer->setPrefix ( multiscaleProblemPrefix + "_Model_" + number2string ( M_ID ) +  label + "_" + number2string ( multiscaleProblemStep - 1 ) );
    importer->setPostDir ( multiscaleProblemFolder );
}

void
MultiscaleModelFSI3D::setExporterFluid ( const IOFilePtr_Type& exporter )
{
    M_fluidVelocity.reset (     new vector_Type ( M_FSIoperator->uFESpace().map(),  M_exporterFluid->mapType() ) );
    M_fluidPressure.reset (     new vector_Type ( M_FSIoperator->pFESpace().map(),  M_exporterFluid->mapType() ) );
    M_fluidDisplacement.reset ( new vector_Type ( M_FSIoperator->mmFESpace().map(), M_exporterFluid->mapType() ) );

    exporter->setMeshProcId ( M_FSIoperator->uFESpace().mesh(), M_FSIoperator->uFESpace().map().comm().MyPID() );

    exporter->addVariable ( IOData_Type::VectorField, "Velocity (fluid)",     M_FSIoperator->uFESpacePtr(),  M_fluidVelocity,     static_cast<UInt> ( 0 ) );
    exporter->addVariable ( IOData_Type::ScalarField, "Pressure (fluid)",     M_FSIoperator->pFESpacePtr(),  M_fluidPressure,     static_cast<UInt> ( 0 ) );
    exporter->addVariable ( IOData_Type::VectorField, "Displacement (fluid)", M_FSIoperator->mmFESpacePtr(), M_fluidDisplacement, static_cast<UInt> ( 0 ) );
}

void
MultiscaleModelFSI3D::setExporterSolid ( const IOFilePtr_Type& exporter )
{
    M_solidDisplacement.reset ( new vector_Type ( M_FSIoperator->dFESpace().map(), M_exporterSolid->mapType() ) );
    M_solidVelocity.reset    ( new vector_Type ( M_FSIoperator->dFESpace().map(), M_exporterSolid->mapType() ) );

    exporter->setMeshProcId ( M_FSIoperator->dFESpace().mesh(), M_FSIoperator->dFESpace().map().comm().MyPID() );

    exporter->addVariable ( IOData_Type::VectorField, "Displacement (solid)", M_FSIoperator->dFESpacePtr(), M_solidDisplacement, static_cast<UInt> ( 0 ) );
    exporter->addVariable ( IOData_Type::VectorField, "Velocity (solid)",     M_FSIoperator->dFESpacePtr(), M_solidVelocity,     static_cast<UInt> ( 0 ) );
}

void
MultiscaleModelFSI3D::exportFluidSolution()
{
    if ( M_FSIoperator->isFluid() )
    {
        M_FSIoperator->exportFluidVelocity ( *M_fluidVelocity );
        M_FSIoperator->exportFluidPressure ( *M_fluidPressure );
        M_FSIoperator->exportFluidDisplacement ( *M_fluidDisplacement );

#ifndef FSI_WITH_EXTERNALPRESSURE
        // Add the external pressure
        *M_fluidPressure += M_externalPressureScalar;
#endif
    }
}

void
MultiscaleModelFSI3D::exportSolidSolution()
{
    if ( M_FSIoperator->isSolid() )
    {
        M_FSIoperator->exportSolidDisplacement ( *M_solidDisplacement );
        M_FSIoperator->exportSolidVelocity ( *M_solidVelocity );
    }
}

void
MultiscaleModelFSI3D::setupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::setupLinearModel() \n";
#endif

    // The linear BCHandler is a copy of the original BCHandler with all BCFunctions giving zero
    M_linearFluidBC.reset ( new bc_Type ( *M_fluidBC->handler() ) );
    M_linearSolidBC.reset ( new bc_Type ( *M_solidBC->handler() ) );

    // Set all the BCFunctions to zero
    BCFunctionBase bcBaseDeltaZero;
    bcBaseDeltaZero.setFunction ( boost::bind ( &MultiscaleModelFSI3D::bcFunctionDeltaZero, this, _1, _2, _3, _4, _5 ) );

    for ( bc_Type::bcBaseIterator_Type i = M_linearFluidBC->begin() ; i != M_linearFluidBC->end() ; ++i )
    {
        i->setBCFunction ( bcBaseDeltaZero );
    }

    for ( bc_Type::bcBaseIterator_Type i = M_linearSolidBC->begin() ; i != M_linearSolidBC->end() ; ++i )
    {
        i->setBCFunction ( bcBaseDeltaZero );
    }

    // Setup linear solution & the RHS
    M_linearSolution.reset ( new vector_Type ( *M_FSIoperator->couplingVariableMap() ) );
    M_linearRHS.reset ( new vector_Type ( *M_FSIoperator->couplingVariableMap() ) );
}

void
MultiscaleModelFSI3D::updateLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::updateLinearModel() \n";
#endif

    //Create the RHS
    *M_linearRHS *= 0;
    M_FSIoperator->bcManageVectorRHS ( M_linearFluidBC, M_linearSolidBC, *M_linearRHS );
}

void
MultiscaleModelFSI3D::solveLinearModel ( bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::solveLinearModel() \n";
#endif

    if ( !solveLinearSystem )
    {
        return;
    }

    imposePerturbation();

    updateLinearModel();

    //Solve the linear problem
    displayModelStatus ( "Solve linear" );
    M_FSIoperator->solveJac ( *M_linearSolution, *M_linearRHS, 0. );

    resetPerturbation();

    //This flag avoid recomputation of the same system
    solveLinearSystem = false;
}

void
MultiscaleModelFSI3D::imposePerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::imposePerturbation() \n";
#endif

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            BCFunctionBase bcBaseDeltaOne;

            multiscaleID_Type boundaryID ( ( *i )->boundaryID ( ( *i )->modelGlobalToLocalID ( M_ID ) ) );
            if ( M_boundaryFlagsAreaPerturbed[boundaryID] == true )
            {
                for ( boundaryAreaFunctionsContainerIterator_Type j = M_boundaryAreaFunctions.begin(); j < M_boundaryAreaFunctions.end(); ++j )
                    if ( ( *j )->fluidFlag() == boundaryFlag ( boundaryID ) )
                    {
                        bcBaseDeltaOne.setFunction ( boost::bind ( &FSI3DBoundaryAreaFunction::functionLinear, *j, _1, _2, _3, _4, _5 ) );
                        M_linearSolidBC->findBCWithFlag ( M_boundaryFlagsArea[boundaryID] ).setBCFunction ( bcBaseDeltaOne );
                    }
            }
            else
            {
                bcBaseDeltaOne.setFunction ( boost::bind ( &MultiscaleModelFSI3D::bcFunctionDeltaOne, this, _1, _2, _3, _4, _5 ) );
                M_linearFluidBC->findBCWithFlag ( boundaryFlag ( boundaryID ) ).setBCFunction ( bcBaseDeltaOne );
            }

            break;
        }
}

void
MultiscaleModelFSI3D::resetPerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8140 ) << "MultiscaleModelFSI3D::resetPerturbation() \n";
#endif

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->isPerturbed() )
        {
            BCFunctionBase bcBaseDeltaZero;
            bcBaseDeltaZero.setFunction ( boost::bind ( &MultiscaleModelFSI3D::bcFunctionDeltaZero, this, _1, _2, _3, _4, _5 ) );

            multiscaleID_Type boundaryID ( ( *i )->boundaryID ( ( *i )->modelGlobalToLocalID ( M_ID ) ) );
            if ( M_boundaryFlagsAreaPerturbed[boundaryID] == true )
            {
                M_linearSolidBC->findBCWithFlag ( M_boundaryFlagsArea[boundaryID] ).setBCFunction ( bcBaseDeltaZero );
            }
            else
            {
                M_linearFluidBC->findBCWithFlag ( boundaryFlag ( boundaryID ) ).setBCFunction ( bcBaseDeltaZero );
                M_boundaryFlagsAreaPerturbed[boundaryID] = true;
            }

            break;
        }
}

} // Namespace multiscale
} // Namespace LifeV
