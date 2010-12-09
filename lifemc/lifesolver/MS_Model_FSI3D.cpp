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

#include <lifemc/lifesolver/MS_Model_FSI3D.hpp>

#include <lifemc/lifesolver/BlockMatrix.hpp>
#include <lifemc/lifesolver/BlockMatrixRN.hpp>
#include <lifemc/lifesolver/ComposedDN.hpp>
#include <lifemc/lifesolver/ComposedNN.hpp>
#include <lifemc/lifesolver/ComposedDNND.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Model_FSI3D::MS_Model_FSI3D() :
        super                          (),
        M_FSIoperator                  ( ),
        M_data                         ( new data_Type() ),
        M_exporterFluid                (),
        M_exporterSolid                (),
        M_importerFluid                (),
        M_importerSolid                (),
        M_solution_tn                  (),
        M_fluidVelocityPressure        (),
        M_fluidDisplacement            (),
        M_solidVelocity                (),
        M_solidDisplacement            (),
        M_fluidBC                      ( new BCInterface_Type() ),
        M_solidBC                      ( new BCInterface_Type() ),
        M_harmonicExtensionBC          ( new BCInterface_Type() ),
//    M_linearizedFluidBC            (),
//    M_linearizedSolidBC            (),
        M_linearBC                     (),
        M_linearRHS                    (),
        M_linearSolution               (),
        M_BCBaseDelta_Zero             (),
        M_BCBaseDelta_One              (),
        M_iter                         (0),
        M_meshDisp_tn                  ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MS_Model_FSI3D::MS_Model_FSI3D() \n";
#endif

    M_type = FSI3D;

    BlockPrecFactory::instance().registerProduct("ComposedDNND",      &ComposedDNND::createComposedDNND);
    BlockPrecFactory::instance().registerProduct("AdditiveSchwarz",   &BlockMatrix::createAdditiveSchwarz) ;
    BlockPrecFactory::instance().registerProduct("AdditiveSchwarzRN", &BlockMatrixRN::createAdditiveSchwarzRN ) ;
    BlockPrecFactory::instance().registerProduct("ComposedDN",        &ComposedDN::createComposedDN ) ;
    BlockPrecFactory::instance().registerProduct("ComposedDN2",       &ComposedDN::createComposedDN2 );

    BlockMatrix::Factory::instance().registerProduct("AdditiveSchwarz",   &BlockMatrix::createAdditiveSchwarz ) ;
    BlockMatrix::Factory::instance().registerProduct("AdditiveSchwarzRN", &BlockMatrixRN::createAdditiveSchwarzRN ) ;

    FSIOperator_Type::solid_raw_type::StructureSolverFactory::instance().registerProduct( "linearVenantKirchhof", &FSIOperator_Type::createLinearStructure );
    //FSIOperator_Type::solid_raw_type::StructureSolverFactory::instance().registerProduct( "nonLinearVenantKirchhof", &FSIOperator_Type::createNonLinearStructure );
}

// ===================================================
// MultiScale PhysicalModel Virtual Methods
// ===================================================
void
MS_Model_FSI3D::SetupData( const std::string& fileName )
{
    super::SetupData( fileName );

    GetPot dataFile( fileName );

    // Load data
    M_data->setup( dataFile );

    if ( M_globalData.get() )
        SetupGlobalData( fileName );

    // Create FSIOperator
    M_FSIoperator = FSIOperator_PtrType( FSIOperator::FSIFactory::instance().createObject( M_data->method() ) );

    // Setup Communicator
    setupCommunicator();

    // Set data
    M_FSIoperator->setData( M_data );
    M_FSIoperator->setDataFile( dataFile );

    // Setup Boundary Conditions for FSI from file
    setupBC( fileName );

    // Setup Boundary Conditions for segregated FSI  from file
//     if( !M_data->isMonolithic() )
//         setupSegregatedBC( fileName );

    // Setup exporter
    if ( M_FSIoperator->isFluid() )
    {
        SetupExporter( M_exporterFluid, dataFile );
        SetupImporter( M_importerFluid, dataFile );
    }
    if ( M_FSIoperator->isSolid() )
    {
        SetupExporter( M_exporterSolid, dataFile, "_Solid" );
        SetupImporter( M_importerSolid, dataFile, "_Solid" );
    }
}

void
MS_Model_FSI3D::SetupModel()
{
    // Mesh transformation (before partitioning, ideally should be done after for scalability)
    //M_FSIoperator->fluidMesh().transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );
    M_FSIoperator->solidMesh().transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    // Mesh partitioning
    M_FSIoperator->partitionMeshes();

    // Mesh transformation (after partitioning - not working for solid)
    M_FSIoperator->fluidMeshPart().mesh()->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );
    //M_FSIoperator->solidMeshPart().mesh()->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );

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
        SetExporterFluid( M_exporterFluid );

    if ( M_FSIoperator->isSolid() )
        SetExporterSolid( M_exporterSolid );

    //Setup solution
    initializeSolution();

    //Setup linear model
    SetupLinearModel();

    M_meshDisp_tn.reset( new vector_Type( *M_fluidDisplacement ));
    //to kill when everithing is merged
    //M_offset=(dynamic_cast<LifeV::Monolithic*>(M_FSIoperator.get()))->getOffset();
}

void
MS_Model_FSI3D::BuildSystem()
{
    // Update BCInterface Operator BC
    updateBC();

    M_FSIoperator->setupSystem();
    M_FSIoperator->buildSystem();
    M_FSIoperator->updateSystem();

    // Update solution at time n
    M_solution_tn.reset( new vector_Type( M_FSIoperator->getSolution() ) );
    M_rhs_tn.reset( new vector_Type( M_FSIoperator->getSolution().getMap() ) );
    M_meshDisp_tn.reset( new vector_Type( M_FSIoperator->meshDisp() ) );
    M_meshDispOld_tn.reset( new vector_Type( M_FSIoperator->meshDisp() ) );

    // Non-linear Richardson solver
    //vector_PtrType solution = M_solution_tn ;
}

void
MS_Model_FSI3D::UpdateSystem()
{
    // Update BCInterface Operator BC
    updateBC();

    // Update System
    M_FSIoperator->updateSystem();


    //Update solution at time n
    *M_meshDispOld_tn = M_FSIoperator->meshMotion().dispOld();
    *M_solution_tn = *M_FSIoperator->solutionPtr();
    *M_meshDisp_tn = M_FSIoperator->meshDisp();
    *M_rhs_tn *= 0;

    boost::dynamic_pointer_cast<Monolithic>(M_FSIoperator)->precPtr()->setRecompute( 1, true );
    //M_solution_tn = M_FSIoperator->solutionPtr();

    M_iter = 0;
}

void
MS_Model_FSI3D::SolveSystem( )
{
    UInt maxSubIterationNumber = M_data->maxSubIterationNumber();
    std::ofstream outRes; // Unuseful variable

    // Non-linear Richardson solver
    //vector_PtrType solution = M_solution_tn ;
    vector_PtrType solution( new vector_Type( *M_solution_tn ) ) ;
    boost::dynamic_pointer_cast<LifeV::Monolithic>(M_FSIoperator)->initializeMesh(M_meshDisp_tn);
    M_FSIoperator->fluid().initialize( *M_fluidVelocityPressure );//useless?
    vector_PtrType newDisp(new vector_Type(*M_solidDisplacement));
    M_FSIoperator->solid().initialize( newDisp, M_solidVelocity );
    vector_PtrType newRHS(new vector_Type(*M_rhs_tn));
    M_FSIoperator->setRHS( newRHS );
    M_FSIoperator->setSolution( *M_solution_tn );
    M_FSIoperator->initializeBDF( *M_solution_tn );
    M_FSIoperator->setRestarts( true );

//    std::cout<< "norm displacement: " << M_meshDisp_tn->Norm2() << std::endl;


    if (M_iter != 0)
    {
        boost::dynamic_pointer_cast<Monolithic>(M_FSIoperator)->precPtr()->setRecompute( 1, false );
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
                                    M_data->dataFluid()->dataTime()->getTime(),
                                    M_iter
                                  );
    if (M_iter == 0)
    {
        *M_fluidDisplacement = M_FSIoperator->meshDisp();
    }


    M_iter=1;
    // We set the solution for the next time step
    M_FSIoperator->setSolutionPtr( solution );

    if ( status == EXIT_FAILURE )
        std::cout << "Non-Linear Richardson failed to converge" << std::endl;

}

void
MS_Model_FSI3D::SaveSolution()
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
        M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->getTime() - M_data->dataFluid()->dataTime()->getTimeStep() );
    if ( M_FSIoperator->isSolid() )
        M_exporterSolid->postProcess( M_data->dataSolid()->dataTime()->getTime() - M_data->dataSolid()->dataTime()->getTimeStep() );

#ifdef HAVE_HDF5
    if ( M_data->dataFluid()->dataTime()->isLastTimeStep() )
    {
        if ( M_FSIoperator->isFluid() )
            ( MS_DynamicCast< hdf5IOFile_Type >( M_exporterFluid ) )->CloseFile();
        if ( M_FSIoperator->isSolid() )
            ( MS_DynamicCast< hdf5IOFile_Type >( M_exporterSolid ) )->CloseFile();
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
MS_Model_FSI3D::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

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
MS_Model_FSI3D::SetupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MS_Model_Fluid3D::SetupLinearModel( ) \n";
#endif

    // Define BCFunctions for tangent problem
    M_BCBaseDelta_Zero.setFunction( boost::bind( &MS_Model_FSI3D::BCFunctionDelta_Zero, this, _1, _2, _3, _4, _5 ) );
    M_BCBaseDelta_One.setFunction(  boost::bind( &MS_Model_FSI3D::BCFunctionDelta_One,  this, _1, _2, _3, _4, _5 ) );

    // The linear BCHandler is a copy of the original BCHandler with all BCFunctions giving zero
    BC_PtrType LinearBCHandler ( new BC_Type( *M_fluidBC->handler() ) );
    M_linearBC = LinearBCHandler;

    // Set all the BCFunctions to zero
    for ( BC_Type::BCBase_Iterator i = M_linearBC->begin() ; i != M_linearBC->end() ; ++i )
        i->setBCFunction( M_BCBaseDelta_Zero );

    // Setup linear solution & the RHS
    M_linearSolution.reset( new vector_Type( M_FSIoperator->un()->getMap() ) );
    M_linearRHS.reset( new vector_Type( M_FSIoperator->un()->getMap() ) );
}

void
MS_Model_FSI3D::UpdateLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MS_Model_FSI3D::UpdateLinearModel() \n";
#endif

    //Create the RHS
    *M_linearRHS *= 0;
    M_FSIoperator->bcManageVectorRHS( M_linearBC, *M_linearRHS );
}

void
MS_Model_FSI3D::SolveLinearModel( bool& SolveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MS_Model_FSI3D::SolveLinearModel() \n";
#endif

    if ( !SolveLinearSystem )
        return;

    ImposePerturbation();

    UpdateLinearModel();

    //Solve the linear problem
    M_FSIoperator->solveJac( *M_linearSolution, *M_linearRHS, 0. );

    ResetPerturbation();

    //This flag avoid recomputation of the same system
    SolveLinearSystem = false;
}

// ===================================================
// Get Methods
// ===================================================
MS_Model_FSI3D::BCInterface_Type&
MS_Model_FSI3D::GetBCInterface()
{
    return *M_fluidBC;
}

Real
MS_Model_FSI3D::GetBoundaryDensity( const BCFlag& /*Flag*/ ) const
{
    return M_FSIoperator->dataFluid()->density();
}

Real
MS_Model_FSI3D::GetBoundaryViscosity( const BCFlag& /*Flag*/ ) const
{
    return M_FSIoperator->dataFluid()->viscosity();
}

Real
MS_Model_FSI3D::GetBoundaryArea( const BCFlag& Flag ) const
{
    return M_FSIoperator->fluid().area( Flag );
}

Real
MS_Model_FSI3D::GetBoundaryFlowRate( const BCFlag& Flag ) const
{
    return M_FSIoperator->fluid().flux( Flag, *M_FSIoperator->solutionPtr() );
}

Real
MS_Model_FSI3D::GetBoundaryPressure( const BCFlag& Flag ) const
{
    return M_FSIoperator->fluid().pressure( Flag, *M_FSIoperator->solutionPtr() );
}

Real
MS_Model_FSI3D::GetBoundaryDynamicPressure( const BCFlag& Flag ) const
{
    return 0.5 * GetBoundaryDensity( Flag ) * ( GetBoundaryFlowRate( Flag ) * GetBoundaryFlowRate( Flag ) )
           / ( GetBoundaryArea( Flag ) * GetBoundaryArea( Flag ) );
}

Real
MS_Model_FSI3D::GetBoundaryLagrangeMultiplier( const BCFlag& Flag ) const
{
    return M_FSIoperator->fluid().LagrangeMultiplier(Flag, *M_fluidBC->handler(), M_FSIoperator->getSolution() );
}

Real
MS_Model_FSI3D::GetBoundaryStress( const BCFlag& Flag, const stressTypes& StressType ) const
{
    switch ( StressType )
    {
    case StaticPressure:
    {
        return -GetBoundaryPressure( Flag );
    }

    case TotalPressure:
    {
        return -GetBoundaryPressure( Flag ) + GetBoundaryDynamicPressure( Flag ) * ( ( GetBoundaryFlowRate( Flag ) > 0.0 ) ? 1 : -1 );
    }

    case LagrangeMultiplier:
    {
        return -GetBoundaryLagrangeMultiplier( Flag );
    }

    default:

        std::cout << "ERROR: Invalid stress type [" << Enum2String( StressType, MS_stressesMap ) << "]" << std::endl;

        return 0.0;
    }
}

Real
MS_Model_FSI3D::GetBoundaryDeltaFlowRate( const BCFlag& Flag, bool& SolveLinearSystem )
{
    SolveLinearModel( SolveLinearSystem );

    return M_FSIoperator->fluid().flux( Flag, *M_linearSolution );
}

Real
MS_Model_FSI3D::GetBoundaryDeltaPressure( const BCFlag& Flag, bool& SolveLinearSystem )
{
    SolveLinearModel( SolveLinearSystem );

    return M_FSIoperator->fluid().pressure( Flag, *M_linearSolution );
}

Real
MS_Model_FSI3D::GetBoundaryDeltaDynamicPressure( const BCFlag& Flag, bool& SolveLinearSystem )
{
    return GetBoundaryDensity( Flag ) * GetBoundaryDeltaFlowRate( Flag, SolveLinearSystem ) * GetBoundaryFlowRate( Flag ) / ( GetBoundaryArea( Flag ) * GetBoundaryArea( Flag ) );
}

Real
MS_Model_FSI3D::GetBoundaryDeltaLagrangeMultiplier( const BCFlag& Flag, bool& SolveLinearSystem )
{
    SolveLinearModel( SolveLinearSystem );

    return M_FSIoperator->fluid().LagrangeMultiplier( Flag, *M_linearBC, *M_linearSolution );
}

Real
MS_Model_FSI3D::GetBoundaryDeltaStress( const BCFlag& Flag, bool& SolveLinearSystem, const stressTypes& StressType )
{
    switch ( StressType )
    {
    case StaticPressure:
    {
        return -GetBoundaryDeltaPressure( Flag, SolveLinearSystem );
    }

    case TotalPressure:
    {
        return -GetBoundaryDeltaPressure( Flag, SolveLinearSystem ) + GetBoundaryDeltaDynamicPressure( Flag, SolveLinearSystem ); //Verify the sign of DynamicPressure contribute!
    }

    case LagrangeMultiplier:
    {
        return -GetBoundaryDeltaLagrangeMultiplier( Flag, SolveLinearSystem );
    }

    default:

        std::cout << "ERROR: Invalid stress type [" << Enum2String( StressType, MS_stressesMap ) << "]" << std::endl;

        return 0.0;
    }
}

// ===================================================
// Private Methods
// ===================================================
void
MS_Model_FSI3D::SetupGlobalData( const std::string& FileName )
{
    GetPot DataFile( FileName );

    //Global data time
    M_data->dataFluid()->setDataTime( M_globalData->GetDataTime() );
    M_data->dataSolid()->setDataTime( M_globalData->GetDataTime() );

    //Global physical quantities
    if ( !DataFile.checkVariable( "fluid/physics/density" ) )
        M_data->dataFluid()->density( M_globalData->GetFluidDensity() );
    if ( !DataFile.checkVariable( "fluid/physics/viscosity" ) )
        M_data->dataFluid()->viscosity( M_globalData->GetFluidViscosity() );

    if ( !DataFile.checkVariable( "solid/physics/density" ) )
        M_data->dataSolid()->setDensity( M_globalData->GetStructureDensity() );
    if ( !DataFile.checkVariable( "solid/physics/poisson" ) )
        M_data->dataSolid()->setPoisson( M_globalData->GetStructurePoissonCoefficient() );
    if ( !DataFile.checkVariable( "solid/physics/young" ) )
        M_data->dataSolid()->setYoung( M_globalData->GetStructureYoungModulus() );

    //M_data->showMe();
}

void
MS_Model_FSI3D::setupCommunicator()
{
    bool fluid = false;
    bool solid = false;

    Int  fluidLeader( 0 );
    Int  solidLeader( 0 );

    boost::shared_ptr<Epetra_Comm> epetraComm;

    if ( M_data->isMonolithic() )
    {
        fluid = true;
        solid = true;
        solidLeader = 0;
        fluidLeader = solidLeader;

        epetraComm = M_comm;
    }
    else
    {
        /*
                MPI_Group  originGroup, newGroup;
                MPI_Comm   newComm;
                MPI_Comm_group( M_comm->Comm(), &originGroup );

                if ( M_comm->NumProc() == 1 )
                {
                    std::cout << "Serial Fluid/Structure computation" << std::endl;
                    newComm = M_comm->Comm();
                    fluid = true;
                    solid = true;
                    fluidLeader = 0;
                    solidLeader = 0;

                    epetraComm = M_comm;
                }
                else
                {
                    Int members[M_comm->NumProc()];

                    solidLeader = 0;
                    fluidLeader = 1-solidLeader;

                    if ( M_comm->MyPID() == solidLeader )
                    {
                        members[0] = solidLeader;
                        MPI_Group_incl( originGroup, 1, members, &newGroup );
                        solid = true;
                    }
                    else
                    {
                        for (Int ii = 0; ii <= M_comm->NumProc(); ++ii)
                        {
                            if ( ii < solidLeader )
                                members[ii] = ii;
                            else if ( ii > solidLeader )
                                members[ii - 1] = ii;
                        }

                        MPI_Group_incl( originGroup, M_comm->NumProc() - 1, members, &newGroup );
                        fluid = true;
                    }

                    MPI_Comm* localComm = new MPI_Comm;
                    MPI_Comm_create( M_comm->Comm(), newGroup, localComm );
                    epetraComm.reset( new Epetra_MpiComm( *localComm ) );
                }
        */
    }

    M_comm->Barrier();

    M_FSIoperator->setFluid( fluid );
    M_FSIoperator->setSolid( solid );

    M_FSIoperator->setFluidLeader( fluidLeader );
    M_FSIoperator->setSolidLeader( solidLeader );

    M_FSIoperator->setComm( epetraComm, M_comm );
}

void
MS_Model_FSI3D::setupBC( const std::string& fileName )
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

/*
void
MS_Model_FSI3D::setupSegregatedBC( const std::string& fileName )
{
    // Setup additional BCHandler for segregated FSI
     M_linearizedFluidBC.reset( new BCInterface_Type() );
     M_linearizedSolidBC.reset( new BCInterface_Type() );
     M_lineardBC.reset( new BCInterface_Type() );

     M_linearizedFluidBC->setPhysicalSolver( M_FSIoperator );
     M_linearizedSolidBC->setPhysicalSolver( M_FSIoperator );

     if ( M_FSIoperator->isFluid() )
         M_linearBC->fillHandler( fileName, "lin_fluid" );

     if ( M_FSIoperator->isFluid() )
         M_linearizedFluidBC->fillHandler( fileName, "lin_fluid" );

     if ( M_FSIoperator->isSolid() )
         M_linearizedSolidBC->fillHandler( fileName, "lin_solid" );

     M_solver->setLinFluidBC( M_linearizedFluidBC->handler() );
     M_solver->setLinSolidBC( M_linearizedSolidBC->handler() );
}
*/

void
MS_Model_FSI3D::updateBC()
{
    if ( M_FSIoperator->isFluid() )
    {
        M_fluidBC->updatePhysicalSolverVariables();
        M_harmonicExtensionBC->updatePhysicalSolverVariables();

//        if( !M_data->isMonolithic() )
//            M_linearBC->updateOperatorVariables();
    }
    else
    {
        M_solidBC->updatePhysicalSolverVariables();

//        if( !M_data->isMonolithic() )
//            M_linearizedSolidBC->updateOperatorVariables();
    }
}

void
MS_Model_FSI3D::SetupExporter( IOFile_PtrType& exporter, const GetPot& dataFile, const std::string& label )
{
    const std::string exporterType = dataFile( "exporter/type", "ensight" );
#ifdef HAVE_HDF5
    if ( exporterType.compare( "hdf5" ) == 0 )
        exporter.reset( new hdf5IOFile_Type( dataFile, label ) );
    else
#endif
        exporter.reset( new ensightIOFile_Type( dataFile, label ) );

    exporter->setPrefix( "Step_" + number2string( MS_ProblemStep ) + "_Model_" + number2string( M_ID ) + label );
    exporter->setDirectory( MS_ProblemFolder );
}

void
MS_Model_FSI3D::SetupImporter( IOFile_PtrType& importer, const GetPot& dataFile, const std::string& label )
{
    const std::string importerType = dataFile( "exporter/type", "ensight" );
#ifdef HAVE_HDF5
    if ( importerType.compare( "hdf5" ) == 0 )
        importer.reset( new hdf5IOFile_Type( dataFile, label ) );
    else
#endif
        importer.reset( new ensightIOFile_Type( dataFile, label ) );

    importer->setPrefix( "Step_" + number2string( MS_ProblemStep - 1 ) + "_Model_" + number2string( M_ID ) + label );
    importer->setDirectory( MS_ProblemFolder );
}

void
MS_Model_FSI3D::SetExporterFluid( const IOFile_PtrType& exporter )
{
    M_fluidVelocityPressure.reset( new vector_Type( M_FSIoperator->fluid().getMap(),  M_exporterFluid->mapType() ) );
    M_fluidDisplacement.reset    ( new vector_Type( M_FSIoperator->mmFESpace().map(), M_exporterFluid->mapType() ) );

    exporter->setMeshProcId( M_FSIoperator->uFESpace().mesh(),
                             M_FSIoperator->uFESpace().map().Comm().MyPID() );

    exporter->addVariable( ExporterData::Vector, "Fluid Velocity", M_fluidVelocityPressure, UInt(0), M_FSIoperator->uFESpace().dof().numTotalDof()  );
    exporter->addVariable( ExporterData::Scalar, "Fluid Pressure", M_fluidVelocityPressure, UInt(3 * M_FSIoperator->uFESpace().dof().numTotalDof()  ),
                           UInt(    M_FSIoperator->pFESpace().dof().numTotalDof()  ) );
    exporter->addVariable( ExporterData::Vector, "Fluid Displacement", M_fluidDisplacement, UInt(0), M_FSIoperator->mmFESpace().dof().numTotalDof() );
}

void
MS_Model_FSI3D::SetExporterSolid( const IOFile_PtrType& exporter )
{
    M_solidDisplacement.reset( new vector_Type( M_FSIoperator->dFESpace().map(), M_exporterSolid->mapType() ) );
    M_solidVelocity.reset    ( new vector_Type( M_FSIoperator->dFESpace().map(), M_exporterSolid->mapType() ) );

    exporter->setMeshProcId( M_FSIoperator->dFESpace().mesh(),
                             M_FSIoperator->dFESpace().map().Comm().MyPID() );

    exporter->addVariable( ExporterData::Vector, "Solid Velocity",     M_solidVelocity,     UInt(0), M_FSIoperator->dFESpace().dof().numTotalDof() );
    exporter->addVariable( ExporterData::Vector, "Solid Displacement", M_solidDisplacement, UInt(0), M_FSIoperator->dFESpace().dof().numTotalDof() );
}

void
MS_Model_FSI3D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MS_Model_FSI3D::InitializeSolution() \n";
#endif

    if ( MS_ProblemStep > 0 )
    {
        M_importerFluid->setMeshProcId( M_FSIoperator->uFESpace().mesh(), M_FSIoperator->uFESpace().map().Comm().MyPID() );
        M_importerSolid->setMeshProcId( M_FSIoperator->dFESpace().mesh(), M_FSIoperator->dFESpace().map().Comm().MyPID() );

        M_importerFluid->addVariable( ExporterData::Vector, "Fluid Velocity", M_fluidVelocityPressure, UInt(0), M_FSIoperator->uFESpace().dof().numTotalDof()  );
        M_importerFluid->addVariable( ExporterData::Scalar, "Fluid Pressure", M_fluidVelocityPressure, UInt(3 * M_FSIoperator->uFESpace().dof().numTotalDof()  ),
                                      UInt(    M_FSIoperator->pFESpace().dof().numTotalDof()  ) );
        M_importerFluid->addVariable( ExporterData::Vector, "Fluid Displacement", M_fluidDisplacement, UInt(0), M_FSIoperator->mmFESpace().dof().numTotalDof() );

        M_importerSolid->addVariable( ExporterData::Vector, "Solid Velocity",     M_solidVelocity,     UInt(0), M_FSIoperator->dFESpace().dof().numTotalDof() );
        M_importerSolid->addVariable( ExporterData::Vector, "Solid Displacement", M_solidDisplacement, UInt(0), M_FSIoperator->dFESpace().dof().numTotalDof() );

        // Import
        M_exporterFluid->setStartIndex( M_importerFluid->importFromTime( M_data->dataFluid()->dataTime()->getInitialTime() ) + 1 );
        M_exporterSolid->setStartIndex( M_importerSolid->importFromTime( M_data->dataSolid()->dataTime()->getInitialTime() ) + 1 );

        // Set Initial solution
        // IMPORTANT NOTE:
        // 1) TODO Remove nRestart flag from the file
        // 2) TODO This part should be rewritten better
        vector_PtrType UniqueVFDOld( new vector_Type( *M_fluidDisplacement, Unique, Zero ) );
        dynamic_cast< Monolithic* >( M_FSIoperator.get())->initializeMesh( UniqueVFDOld );

        vector_PtrType initSol( new EpetraVector( *M_FSIoperator->getCouplingVariableMap() ) );
        vector_PtrType UniqueV( new vector_Type( *M_fluidVelocityPressure, Unique, Zero ) );
        *initSol = *UniqueV;
        M_FSIoperator->fluid().initialize( *initSol );

        UniqueV.reset( new vector_Type( *M_FSIoperator->getCouplingVariableMap(), Unique, Zero ) );
        UInt offset = dynamic_cast<Monolithic*>(M_FSIoperator.get())->getOffset();
        UniqueV->subset( *M_solidDisplacement, M_solidDisplacement->getMap(), (UInt)0, offset );
        *UniqueV *= 1 / ( M_FSIoperator->solid().rescaleFactor() * M_data->dataFluid()->dataTime()->getTimeStep() );
        M_FSIoperator->solid().initialize( UniqueV );
        *initSol += *UniqueV;

        if ( !M_data->method().compare("monolithicGI") )
        {
            vector_PtrType UniqueVFD ( new vector_Type( *M_FSIoperator->getCouplingVariableMap(), Unique, Zero ) );
            UniqueVFD->subset( *M_fluidDisplacement, M_fluidDisplacement->getMap(), (UInt)0, dynamic_cast<MonolithicGI*>(M_FSIoperator.get())->mapWithoutMesh().getMap(Unique)->NumGlobalElements());
            *initSol += *UniqueVFD;
        }

        vector_PtrType initSolSVel( new vector_Type( *M_FSIoperator->getCouplingVariableMap(), Unique, Zero ) );
        initSolSVel->subset( *M_solidVelocity,M_solidVelocity->getMap(), (UInt)0, offset );
        *initSolSVel *= 1 / ( M_FSIoperator->solid().rescaleFactor() * M_data->dataSolid()->dataTime()->getTimeStep() );
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
MS_Model_FSI3D::ImposePerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MS_Model_FSID::ImposePerturbation() \n";
#endif

    for ( MS_CouplingsVector_ConstIterator i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->IsPerturbed() )
        {
            M_linearBC->GetBCWithFlag( ( *i )->GetFlag( ( *i )->GetModelLocalID( M_ID ) ) ).setBCFunction( M_BCBaseDelta_One );

            break;
        }
}

void
MS_Model_FSI3D::ResetPerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8140 ) << "MS_Model_FSI3D::ResetPerturbation() \n";
#endif

    for ( MS_CouplingsVector_ConstIterator i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->IsPerturbed() )
        {
            M_linearBC->GetBCWithFlag( ( *i )->GetFlag( ( *i )->GetModelLocalID( M_ID ) ) ).setBCFunction( M_BCBaseDelta_Zero );

            break;
        }
}

Real
MS_Model_FSI3D::BCFunctionDelta_Zero( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/ )
{
    return 0.;
}

Real
MS_Model_FSI3D::BCFunctionDelta_One( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/ )
{
    return 1.;
}

} // Namespace LifeV
