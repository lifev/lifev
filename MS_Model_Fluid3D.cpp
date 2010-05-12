//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief MultiScale Model Fluid3D
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 12-03-2009
 */

#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Model_Fluid3D::MS_Model_Fluid3D() :
    super                          (),
    M_exporter                     ( new IOFile_Type() ),
    M_importer                     ( new IOFile_Type() ),
    M_FileName                     (),
    M_Fluid                        (),
    M_FluidBC                      ( new BCInterface_Type() ),
    M_FluidBDF                     (),
    M_FluidData                    ( new Data_Type() ),
    M_FluidMesh                    (),
    M_FluidFullMap                 (),
    M_FluidSolution                (),
    M_LinearFluidBC                ( new BCInterface_Type() ),
    M_UpdateLinearModel            ( true ),
    M_uFESpace                     (),
    M_pFESpace                     (),
    M_uDOF                         ( 0 ),
    M_pDOF                         ( 0 ),
    M_lmDOF                        ( 0 ),
    M_alpha                        ( 0 ),
    M_beta                         (),
    M_RHS                          (),
    M_SubiterationsMaximumNumber   (),
    M_Tolerance                    (),
    M_generalizedAitken            ()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::MS_Model_Fluid3D() \n";
#endif

    M_type = Fluid3D;
}

MS_Model_Fluid3D::MS_Model_Fluid3D( const MS_Model_Fluid3D& Fluid3D ) :
    super                          ( Fluid3D ),
    M_exporter                     ( Fluid3D.M_exporter ),
    M_importer                     ( Fluid3D.M_importer ),
    M_FileName                     ( Fluid3D.M_FileName ),
    M_Fluid                        ( Fluid3D.M_Fluid ),
    M_FluidBC                      ( Fluid3D.M_FluidBC ),
    M_FluidBDF                     ( Fluid3D.M_FluidBDF ),
    M_FluidData                    ( Fluid3D.M_FluidData ),
    M_FluidMesh                    ( Fluid3D.M_FluidMesh ),
    M_FluidFullMap                 ( Fluid3D.M_FluidFullMap ),
    M_FluidSolution                ( Fluid3D.M_FluidSolution ),
    M_LinearFluidBC                ( Fluid3D.M_LinearFluidBC ),
    M_UpdateLinearModel            ( Fluid3D.M_UpdateLinearModel ),
    M_uFESpace                     ( Fluid3D.M_uFESpace ),
    M_pFESpace                     ( Fluid3D.M_pFESpace ),
    M_uDOF                         ( Fluid3D.M_uDOF ),
    M_pDOF                         ( Fluid3D.M_pDOF ),
    M_lmDOF                        ( Fluid3D.M_lmDOF ),
    M_alpha                        ( Fluid3D.M_alpha ),
    M_beta                         ( Fluid3D.M_beta ),
    M_RHS                          ( Fluid3D.M_RHS ),
    M_SubiterationsMaximumNumber   ( Fluid3D.M_SubiterationsMaximumNumber ),
    M_Tolerance                    ( Fluid3D.M_Tolerance ),
    M_generalizedAitken            ( Fluid3D.M_generalizedAitken )
{
}

// ===================================================
// Operators
// ===================================================
MS_Model_Fluid3D&
MS_Model_Fluid3D::operator=( const MS_Model_Fluid3D& Fluid3D )
{
    if ( this != &Fluid3D )
    {
        super::operator=( Fluid3D );
        M_exporter                     = Fluid3D.M_exporter;
        M_importer                     = Fluid3D.M_importer;
        M_FileName                     = Fluid3D.M_FileName;
        M_Fluid                        = Fluid3D.M_Fluid;
        M_FluidBC                      = Fluid3D.M_FluidBC;
        M_FluidBDF                     = Fluid3D.M_FluidBDF;
        M_FluidMesh                    = Fluid3D.M_FluidMesh;
        M_FluidFullMap                 = Fluid3D.M_FluidFullMap;
        M_FluidSolution                = Fluid3D.M_FluidSolution;
        M_LinearFluidBC                = Fluid3D.M_LinearFluidBC;
        M_UpdateLinearModel            = Fluid3D.M_UpdateLinearModel;
        M_uFESpace                     = Fluid3D.M_uFESpace;
        M_pFESpace                     = Fluid3D.M_pFESpace;
        M_uDOF                         = Fluid3D.M_uDOF;
        M_pDOF                         = Fluid3D.M_pDOF;
        M_lmDOF                        = Fluid3D.M_lmDOF;
        M_alpha                        = Fluid3D.M_alpha;
        M_beta                         = Fluid3D.M_beta;
        M_RHS                          = Fluid3D.M_RHS;
        M_SubiterationsMaximumNumber   = Fluid3D.M_SubiterationsMaximumNumber;
        M_Tolerance                    = Fluid3D.M_Tolerance;
        M_generalizedAitken            = Fluid3D.M_generalizedAitken;
    }

    return *this;
}

// ===================================================
// MultiScale PhysicalModel Virtual Methods
// ===================================================
void
MS_Model_Fluid3D::SetupData( const std::string& FileName )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupData( ) \n";
#endif

    super::SetupData( FileName );
    M_FileName = FileName;

    GetPot DataFile( FileName );

    //Fluid data
    M_FluidData->setDataTime( M_dataPhysics->GetDataTime() );
    M_FluidData->setup( DataFile );

    //If global physical quantities are specified also in the data file, replace them
    if ( !DataFile.checkVariable( "fluid/physics/density" ) )
        M_FluidData->density( M_dataPhysics->GetFluidDensity() );
    if ( !DataFile.checkVariable( "fluid/physics/viscosity" ) )
        M_FluidData->viscosity( M_dataPhysics->GetFluidViscosity() );

    // Parameters for the NS Iterations
    M_SubiterationsMaximumNumber = DataFile( "fluid/miscellaneous/SubITMax", 0 );
    M_Tolerance                  = DataFile( "fluid/miscellaneous/Tolerance", 1.e-6 );

    M_generalizedAitken.setDefault(          DataFile( "fluid/miscellaneous/Omega",         1.e-3 ) );
    M_generalizedAitken.UseDefaultOmega(     DataFile( "fluid/miscellaneous/fixedOmega",    false ) );
    M_generalizedAitken.setMinimizationType( DataFile( "fluid/miscellaneous/inverseOmega", true ) );

    //Boundary Conditions for the problem
    M_FluidBC->SetOperator( M_Fluid );    //MUST BE MOVED AFTER M_Fluid.reset !!!
    M_FluidBC->FillHandler( FileName, "fluid" );
    M_FluidBC->UpdateOperatorVariables(); //MUST BE MOVED INSIDE THE UPDATE !!!

    //Setup linear problem
    SetupLinearData( FileName );

    //Exporter
    M_exporter->setDataFromGetPot( DataFile );
    M_exporter->setPrefix( "Step_" + number2string( MS_ProblemStep ) + "_Model_" + number2string( M_ID ) );
    M_exporter->setDirectory( MS_ProblemFolder );

    //Importer
    M_importer->setDataFromGetPot( DataFile );
    M_importer->setPrefix( "Step_" + number2string( MS_ProblemStep - 1 ) + "_Model_" + number2string( M_ID ) );
    M_importer->setDirectory( MS_ProblemFolder );
}

void
MS_Model_Fluid3D::SetupModel()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupProblem() \n";
#endif

    //FEspace
    SetupFEspace();

    //DOF
    SetupDOF();

    //Add flow rate offset to BC
    SetupBCOffset( M_FluidBC->GetHandler() );

    //Fluid
    M_Fluid.reset( new Fluid_Type( *M_FluidData, *M_uFESpace, *M_pFESpace, *M_comm, M_lmDOF ) );
    GetPot DataFile( M_FileName );
    M_Fluid->setUp( DataFile ); //Remove Preconditioner and Solver if possible!

    //Fluid MAP
    M_FluidFullMap.reset( new EpetraMap( M_Fluid->getMap() ) );

    //BDF
    M_FluidBDF.reset( new BDF_Type( M_FluidData->dataTime()->getBDF_order() ) );

    //Problem coefficients
    M_beta.reset( new VectorType( M_FluidFullMap ) );
    M_RHS.reset ( new VectorType( M_FluidFullMap ) );

    //Post-processing
    M_exporter->setMeshProcId( M_FluidMesh->mesh(), M_comm->MyPID() );

    //M_FluidSolution.reset( new VectorType( M_FluidFullMap, Repeated ) );
    //M_FluidSolution.reset( new VectorType( M_Fluid->solution(), Repeated ) );

    M_FluidSolution.reset( new VectorType( M_Fluid->solution(), M_exporter->mapType() ) );
    if ( M_exporter->mapType() == Unique )
        M_FluidSolution->setCombineMode( Zero );

    M_exporter->addVariable( ExporterData::Vector, "Velocity", M_FluidSolution, static_cast <UInt> ( 0 ), M_uDOF );
#ifdef HAVE_HDF5
    M_exporter->addVariable( ExporterData::Scalar, "Pressure", M_FluidSolution, 3 * M_uDOF, M_pDOF);
#else
    M_exporter->addVariable( ExporterData::Scalar, "Pressure", M_FluidSolution, 3 * M_uDOF, 3 * M_uDOF + M_pDOF );
#endif

    //Setup linear model
    SetupLinearModel();

    //Setup solution
    SetupSolution();
}

void
MS_Model_Fluid3D::BuildSystem()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::BuildSystem() \n";
#endif

    //Build constant matrices
    M_Fluid->buildSystem();

    //Initialize BDF
    M_FluidBDF->bdf_u().initialize_unk( M_Fluid->solution() );

    //Define problem coefficients
    if ( M_FluidData->Stokes() )
    {
        M_alpha  = 0.0;
        *M_beta  = M_Fluid->solution();
        *M_RHS  *= 0.0;
    }
    else
    {
        M_alpha = M_FluidBDF->bdf_u().coeff_der( 0 ) / M_FluidData->dataTime()->getTimeStep();
        *M_beta = M_FluidBDF->bdf_u().extrap();
        *M_RHS  = M_Fluid->matrMass() * M_FluidBDF->bdf_u().time_der( M_FluidData->dataTime()->getTimeStep() );
    }

    //Set problem coefficients
    M_Fluid->updateSystem( M_alpha, *M_beta, *M_RHS );
}

void
MS_Model_Fluid3D::UpdateSystem()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::UpdateSystem() \n";
#endif

    //Update BDF
    M_FluidBDF->bdf_u().shift_right( M_Fluid->solution() );

    //Update problem coefficients
    M_alpha = M_FluidBDF->bdf_u().coeff_der( 0 ) / M_FluidData->dataTime()->getTimeStep();
    *M_beta = M_FluidBDF->bdf_u().extrap();
    *M_RHS  = M_Fluid->matrMass() * M_FluidBDF->bdf_u().time_der( M_FluidData->dataTime()->getTimeStep() );

    //Set problem coefficients
    M_Fluid->updateSystem( M_alpha, *M_beta, *M_RHS );

    //Linear system need to be updated
    M_UpdateLinearModel = true;

    //Recompute preconditioner
    M_Fluid->resetPrec( true );
}

void
MS_Model_Fluid3D::SolveSystem()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SolveSystem() \n";
#endif

    //Solve the problem
    M_Fluid->iterate( *M_FluidBC->GetHandler() );

    if ( !M_FluidData->Stokes() )
    {
        Real residual = ( *M_beta - M_Fluid->solution() ).Norm2(); // Residual is computed on the whole solution vector;

        if ( M_displayer->isLeader() )
            std::cout << "  F-  Residual:                                " << residual << std::endl;

        M_generalizedAitken.restart( true );
        for ( UInt subIT = 1; subIT <= M_SubiterationsMaximumNumber; ++subIT )
        {
            *M_beta += M_generalizedAitken.computeDeltaLambdaScalar( *M_beta, *M_beta - M_Fluid->solution() );

            //Linear model need to be updated!
            M_Fluid->updateSystem( M_alpha, *M_beta, *M_RHS );
            M_UpdateLinearModel = true;

            //Solve system
            M_Fluid->iterate( *M_FluidBC->GetHandler() );

            // Check the new residual
            residual = ( *M_beta - M_Fluid->solution() ).Norm2(); // Residual is computed on the whole solution vector

            // Display subiteration information
            if ( M_displayer->isLeader() )
            {
                std::cout << "  F-  Sub-iteration n.:                        " << subIT << std::endl;
                std::cout << "  F-  Residual:                                " << residual << std::endl;
            }

            // Verify tolerance
            if ( residual <= M_Tolerance )
                break;
        }
    }
}

void
MS_Model_Fluid3D::SaveSolution()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SaveSolution() \n";
#endif

    //Post-processing
    *M_FluidSolution = M_Fluid->solution();
    M_exporter->postProcess( M_FluidData->dataTime()->getTime() );

#ifdef HAVE_HDF5
    if ( M_FluidData->dataTime()->isLastTimeStep() )
        M_exporter->CloseFile();
#endif

}

void
MS_Model_Fluid3D::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "uOrder              = " << M_FluidData->uOrder() << std::endl
                  << "pOrder              = " << M_FluidData->pOrder() << std::endl << std::endl;

        std::cout << "uDOF                = " << 3 * M_uDOF << std::endl
                  << "pDOF                = " << M_pDOF << std::endl
                  << "lmDOF               = " << M_lmDOF << std::endl << std::endl;

        std::cout << "maxH                = " << M_FluidData->dataMesh()->mesh()->maxH() << std::endl
                  << "meanH               = " << M_FluidData->dataMesh()->mesh()->meanH() << std::endl << std::endl;

        std::cout << "NS SubITMax         = " << M_SubiterationsMaximumNumber << std::endl
                  << "NS Tolerance        = " << M_Tolerance << std::endl << std::endl << std::endl << std::endl;
    }
}

// ===================================================
// Methods
// ===================================================
void
MS_Model_Fluid3D::SetupLinearData( const std::string& FileName )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupLinearData( ) \n";
#endif

    // Boundary Conditions for the linear problem
    M_LinearFluidBC->SetOperator( M_Fluid );
    M_LinearFluidBC->FillHandler( FileName, "linear_fluid" );
}

void
MS_Model_Fluid3D::SetupLinearModel()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupLinearModel( ) \n";
#endif

    SetupBCOffset( M_LinearFluidBC->GetHandler() );
}

void
MS_Model_Fluid3D::UpdateLinearModel()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::UpdateLinearModel() \n";
#endif

    //Create an empty vector
    VectorType VectorZero( *M_FluidSolution ); VectorZero = 0.0;

    //UpdateLinearModel
    M_Fluid->updateLinearSystem( M_Fluid->matrNoBC(),
                                 M_alpha,
                                 *M_beta,
                                 M_Fluid->solution(),
                                 VectorZero,
                                 VectorZero,
                                 VectorZero,
                                 VectorZero );

    //Update Properties of BC
    M_LinearFluidBC->UpdateOperatorVariables();

    //Linear System Updated
    M_UpdateLinearModel = false;
}

void
MS_Model_Fluid3D::SolveLinearModel( bool& SolveLinearSystem )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SolveLinearModel() \n";
#endif

    if ( !SolveLinearSystem )
        return;

    if ( M_UpdateLinearModel )
        UpdateLinearModel();

    //Solve the linear problem
    M_Fluid->iterateLin( *M_LinearFluidBC->GetHandler() );

    //This flag avoid recomputation of the same system
    SolveLinearSystem = false;
}

// ===================================================
// Get Methods
// ===================================================
MS_Model_Fluid3D::BCInterface_Type&
MS_Model_Fluid3D::GetBCInterface()
{
    return *M_FluidBC;
}

MS_Model_Fluid3D::BCInterface_Type&
MS_Model_Fluid3D::GetLinearBCInterface()
{
    return *M_LinearFluidBC;
}

Real
MS_Model_Fluid3D::GetBoundaryDensity( const BCFlag& /*Flag*/ ) const
{
    return M_FluidData->density();
}

Real
MS_Model_Fluid3D::GetBoundaryViscosity( const BCFlag& /*Flag*/ ) const
{
    return M_FluidData->viscosity();
}

Real
MS_Model_Fluid3D::GetBoundaryArea( const BCFlag& Flag ) const
{
    return M_Fluid->area( Flag );
}

Real
MS_Model_Fluid3D::GetBoundaryFlowRate( const BCFlag& Flag ) const
{
    return M_Fluid->flux( Flag );
}

Real
MS_Model_Fluid3D::GetBoundaryPressure( const BCFlag& Flag ) const
{
    return M_Fluid->pressure( Flag );
}

Real
MS_Model_Fluid3D::GetBoundaryDynamicPressure( const BCFlag& Flag ) const
{
    return 0.5 * GetBoundaryDensity( Flag ) * ( GetBoundaryFlowRate( Flag ) * GetBoundaryFlowRate( Flag ) )
                                            / ( GetBoundaryArea( Flag ) * GetBoundaryArea( Flag ) );
}

Real
MS_Model_Fluid3D::GetBoundaryLagrangeMultiplier( const BCFlag& Flag ) const
{
    return M_Fluid->LagrangeMultiplier(Flag, *M_FluidBC->GetHandler() );
}

Real
MS_Model_Fluid3D::GetBoundaryStress( const BCFlag& Flag, const stressTypes& StressType ) const
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

            std::cout << "ERROR: Invalid stress type [" << Enum2String( StressType, stressMap ) << "]" << std::endl;

            return 0.0;
    }
}

Real
MS_Model_Fluid3D::GetBoundaryDeltaFlux( const BCFlag& Flag, bool& SolveLinearSystem )
{
    SolveLinearModel( SolveLinearSystem );

    return M_Fluid->GetLinearFlux( Flag );
}

Real
MS_Model_Fluid3D::GetBoundaryDeltaPressure( const BCFlag& Flag, bool& SolveLinearSystem )
{
    SolveLinearModel( SolveLinearSystem );

    return M_Fluid->GetLinearPressure( Flag );
}

Real
MS_Model_Fluid3D::GetBoundaryDeltaDynamicPressure( const BCFlag& Flag, bool& SolveLinearSystem )
{
    SolveLinearModel( SolveLinearSystem );

    return GetBoundaryDensity( Flag ) * M_Fluid->GetLinearFlux( Flag ) * GetBoundaryFlowRate( Flag ) / ( GetBoundaryArea( Flag ) * GetBoundaryArea( Flag ) );
}

Real
MS_Model_Fluid3D::GetBoundaryDeltaLagrangeMultiplier( const BCFlag& Flag, bool& SolveLinearSystem )
{
    SolveLinearModel( SolveLinearSystem );

    return M_Fluid->LinearLagrangeMultiplier(Flag, *M_LinearFluidBC->GetHandler() );
}

Real
MS_Model_Fluid3D::GetBoundaryDeltaStress( const BCFlag& Flag, bool& SolveLinearSystem, const stressTypes& StressType )
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

            std::cout << "ERROR: Invalid stress type [" << Enum2String( StressType, stressMap ) << "]" << std::endl;

            return 0.0;
    }
}

// ===================================================
// Private Methods
// ===================================================
void
MS_Model_Fluid3D::SetupFEspace()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupFEspace() \n";
#endif

    //Transform mesh
    M_FluidData->dataMesh()->mesh()->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    //Partition mesh
    M_FluidMesh.reset( new PartitionMesh_Type( *M_FluidData->dataMesh()->mesh(), *M_comm ) );
    M_FluidData->dataMesh()->setMesh( M_FluidMesh->mesh() );

    //Velocity FE Space
    const RefFE* u_refFE;
    const QuadRule* u_qR;
    const QuadRule* u_bdQr;

    if ( M_FluidData->uOrder().compare( "P2" ) == 0 )
    {
        u_refFE = &feTetraP2;
        u_qR = &quadRuleTetra15pt; // DoE 5
        u_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else
        if ( M_FluidData->uOrder().compare( "P1" ) == 0 )
        {
            u_refFE = &feTetraP1;
            u_qR = &quadRuleTetra4pt; // DoE 2
            u_bdQr = &quadRuleTria3pt; // DoE 2
        }
        else
            if ( M_FluidData->uOrder().compare( "P1Bubble" ) == 0 )
            {
                u_refFE = &feTetraP1bubble;
                u_qR = &quadRuleTetra64pt; // DoE 2
                u_bdQr = &quadRuleTria3pt; // DoE 2
            }
            else
            {
                if ( M_displayer->isLeader() )
                    std::cout << M_FluidData->uOrder() << " Velocity FE not implemented yet." << std::endl;
                exit( EXIT_FAILURE );
            }

    //Pressure FE Space
    const RefFE* p_refFE;
    const QuadRule* p_qR;
    const QuadRule* p_bdQr;

    if ( M_FluidData->pOrder().compare( "P2" ) == 0 )
    {
        p_refFE = &feTetraP2;
        p_qR = u_qR;
        p_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else
        if ( M_FluidData->pOrder().compare( "P1" ) == 0 )
        {
            p_refFE = &feTetraP1;
            p_qR = u_qR;
            p_bdQr = &quadRuleTria3pt; // DoE 2
        }
        else
        {
            if ( M_displayer->isLeader() )
                std::cout << M_FluidData->pOrder() << " pressure FE not implemented yet." << std::endl;
            exit( EXIT_FAILURE );
        }

    M_uFESpace.reset( new FESpace_Type( *M_FluidMesh, *u_refFE, *u_qR, *u_bdQr, 3, *M_comm ) );
    M_pFESpace.reset( new FESpace_Type( *M_FluidMesh, *p_refFE, *p_qR, *p_bdQr, 1, *M_comm ) );
}

void
MS_Model_Fluid3D::SetupDOF()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupDOF \n";
#endif

    //DOF
    M_uDOF  = M_uFESpace->dof().numTotalDof();
    M_pDOF  = M_pFESpace->dof().numTotalDof();
    M_lmDOF = M_FluidBC->GetHandler()->getNumberBCWithType( Flux );

    //M_uDOF = M_uFESpace->map().getMap(Unique)->NumGlobalElements();
    //M_pDOF = M_pFESpace->map().getMap(Unique)->NumGlobalElements();
}

void
MS_Model_Fluid3D::SetupBCOffset( const boost::shared_ptr< BC_Type >& BC )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupBCOffset( BC ) \n";
#endif

    UInt offset = M_uFESpace->map().getMap( Unique )->NumGlobalElements() + M_pFESpace->map().getMap( Unique )->NumGlobalElements();

    std::vector< BCName > FluxVector = BC->getBCWithType( Flux );
    for ( UInt i = 0; i < M_lmDOF; ++i )
        BC->setOffset( FluxVector[i], offset + i );
}

void
MS_Model_Fluid3D::SetupSolution()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupSolution() \n";
#endif

    if ( MS_ProblemStep > 0 )
    {
        M_importer->setMeshProcId( M_FluidMesh->mesh(), M_comm->MyPID() );

        M_importer->addVariable( ExporterData::Vector, "Velocity", M_FluidSolution, static_cast <UInt> ( 0 ), M_uDOF );
        #ifdef HAVE_HDF5
            M_importer->addVariable( ExporterData::Scalar, "Pressure", M_FluidSolution, 3 * M_uDOF, M_pDOF);
        #else
            M_importer->addVariable( ExporterData::Scalar, "Pressure", M_FluidSolution, 3 * M_uDOF, 3 * M_uDOF + M_pDOF );
        #endif

        // Import
        M_exporter->setStartIndex( M_importer->importFromTime( M_FluidData->dataTime()->getInitialTime() ) + 1 );
    }
    else
        *M_FluidSolution = 0.0;

    M_Fluid->initialize( *M_FluidSolution );
}

} // Namespace LifeV
