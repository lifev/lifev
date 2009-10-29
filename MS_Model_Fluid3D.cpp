/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-03-12

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/**
 \file MS_Model_Fluid3D.cpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-03-12
 */

#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>

namespace LifeV {

// ===================================================
//! Constructor
// ===================================================
MS_Model_Fluid3D::MS_Model_Fluid3D() :
    super           (),
    M_output        (),
    M_Fluid         (),
    M_FluidBC       (),
    M_FluidBDF      (),
    M_FluidData     (),
    M_FluidMesh     (),
    M_FluidFullMap  (),
    M_FluidSolution (),
    M_LinearFluidBC (),
    M_uFESpace      (),
    M_pFESpace      (),
    M_uDOF          ( 0 ),
    M_pDOF          ( 0 ),
    M_alpha         ( 0 ),
    M_beta          (),
    M_RHS           ()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::MS_Model_Fluid3D() \n";
#endif

    M_type = Fluid3D;
}

MS_Model_Fluid3D::MS_Model_Fluid3D( const MS_Model_Fluid3D& Fluid3D ) :
    super           ( Fluid3D ),
    M_output        ( Fluid3D.M_output ),
    M_Fluid         ( Fluid3D.M_Fluid ),
    M_FluidBC       ( Fluid3D.M_FluidBC ),
    M_FluidBDF      ( Fluid3D.M_FluidBDF ),
    M_FluidData     ( Fluid3D.M_FluidData ),
    M_FluidMesh     ( Fluid3D.M_FluidMesh ),
    M_FluidFullMap  ( Fluid3D.M_FluidFullMap ),
    M_FluidSolution ( Fluid3D.M_FluidSolution ),
    M_LinearFluidBC ( Fluid3D.M_LinearFluidBC ),
    M_uFESpace      ( Fluid3D.M_uFESpace ),
    M_pFESpace      ( Fluid3D.M_pFESpace ),
    M_uDOF          ( Fluid3D.M_uDOF ),
    M_pDOF          ( Fluid3D.M_pDOF ),
    M_alpha         ( Fluid3D.M_alpha ),
    M_beta          ( Fluid3D.M_beta ),
    M_RHS           ( Fluid3D.M_RHS )
{
}

// ===================================================
//! Methods
// ===================================================
MS_Model_Fluid3D&
MS_Model_Fluid3D::operator=( const MS_Model_Fluid3D& Fluid3D )
{
    if ( this != &Fluid3D )
    {
        super::operator=( Fluid3D );
        M_output        = Fluid3D.M_output;
        M_Fluid         = Fluid3D.M_Fluid;
        M_FluidBC       = Fluid3D.M_FluidBC;
        M_FluidBDF      = Fluid3D.M_FluidBDF;
        M_FluidMesh     = Fluid3D.M_FluidMesh;
        M_FluidFullMap  = Fluid3D.M_FluidFullMap;
        M_FluidSolution = Fluid3D.M_FluidSolution;
        M_LinearFluidBC = Fluid3D.M_LinearFluidBC;
        M_uFESpace      = Fluid3D.M_uFESpace;
        M_pFESpace      = Fluid3D.M_pFESpace;
        M_uDOF          = Fluid3D.M_uDOF;
        M_pDOF          = Fluid3D.M_pDOF;
        M_alpha         = Fluid3D.M_alpha;
        M_beta          = Fluid3D.M_beta;
        M_RHS           = Fluid3D.M_RHS;
    }

    return *this;
}

void
MS_Model_Fluid3D::SetupLinearData( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupLinearData( ) \n";
#endif

    // Boundary Conditions for the linear problem
    M_LinearFluidBC.reset( new FluidBCType( M_dataFile, "linear_fluid" ) );
    M_LinearFluidBC->SetOperator( M_Fluid );
    M_LinearFluidBC->BuildHandler();
}

void
MS_Model_Fluid3D::SetupLinearModel( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupLinearModel( ) \n";
#endif

    // Impose offset for Flux BC
    imposeFluxes( M_LinearFluidBC );
}

void
MS_Model_Fluid3D::BuildLinearSystem( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::BuildLinearSystem() \n";
#endif

    UpdateLinearSystem();
}

void
MS_Model_Fluid3D::UpdateLinearSystem( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::UpdateLinearSystem() \n";
#endif

    //Create an empty vector
    FluidVectorType VectorZero( *M_FluidSolution ); VectorZero = 0.0;

    //updateLinearSystem
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

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Model_Fluid3D::SolveLinearSystem( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SolveLinearSystem() \n";
#endif

    //Solve the linear problem
    M_Fluid->iterateLin( *M_LinearFluidBC->Handler_ptr() );

    //MPI Barrier
    M_comm->Barrier();
}

// ===================================================
//! MultiScale Physical Model
// ===================================================
void
MS_Model_Fluid3D::SetupData( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupData( ) \n";
#endif

    //Fluid data
    M_FluidData.reset( new FluidDataType( M_dataFile ) );

    //dataPhysics
    if ( !M_dataFile.checkVariable( "fluid/physics/density" ) )
        M_FluidData->density( M_dataPhysics->GetFluidDensity() );
    if ( !M_dataFile.checkVariable( "fluid/physics/viscosity" ) )
        M_FluidData->viscosity( M_dataPhysics->GetFluidViscosity() );

    //dataTime
    M_FluidData->setInitialTime( M_dataTime->getInitialTime() );
    M_FluidData->setEndTime( M_dataTime->getEndTime() );

    M_FluidData->setTime( M_dataTime->getInitialTime() );
    M_FluidData->setTimeStep( M_dataTime->getTimeStep() );
    //M_FluidData->setTimeStep	( globalTimeStep / std::ceil( globalTimeStep / M_FluidData->getTimeStep() ) );

    //FEspace & DOF
    setupFEspace();
    setupDOF();

    //Boundary Conditions for the problem
    M_FluidBC.reset( new FluidBCType( M_dataFile, "fluid" ) );
    M_FluidBC->SetOperator( M_Fluid ); //MUST BE MOVED AFTER M_Fluid.reset !!!
    M_FluidBC->BuildHandler();
    M_FluidBC->UpdateOperatorVariables(); //MUST BE MOVED INSIDE THE UPDATE !!!

    SetupLinearData(); //TODO MOVE THIS!
}

void
MS_Model_Fluid3D::SetupModel( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupProblem() \n";
#endif

    //Fluid
    M_Fluid.reset( new FluidType( *M_FluidData, *M_uFESpace, *M_pFESpace, *M_comm, imposeFluxes( M_FluidBC ) ) );

    M_Fluid->setUp( M_dataFile ); //Remove Preconditioner and Solver if possible!

    //Fluid MAP
    M_FluidFullMap.reset( new EpetraMap( M_Fluid->getMap() ) );

    //BDF
    M_FluidBDF.reset( new FluidBDFType( M_FluidData->getBDF_order() ) );

    //Problem coefficients
    M_beta.reset( new FluidVectorType( M_FluidFullMap ) );
    M_RHS.reset ( new FluidVectorType( M_FluidFullMap ) );

    //Post-processing
    M_output.reset( new OutputType( M_dataFile, M_FluidMesh->mesh(), M_modelName + "_" + number2string( M_ID ), M_comm->MyPID() ) );

    M_FluidSolution.reset( new FluidVectorType( M_FluidFullMap, Repeated ) );
    //M_FluidSolution.reset( new FluidVectorType( M_Fluid->solution(), Repeated ) );

    M_output->addVariable( ExporterData::Vector, "velocity", M_FluidSolution, UInt( 0 ), M_uDOF );
    M_output->addVariable( ExporterData::Scalar, "pressure", M_FluidSolution, 3 * M_uDOF, 3 * M_uDOF + M_pDOF );

    //MPI Barrier
    M_comm->Barrier();

    SetupLinearModel(); //TODO MOVE THIS!
}

void
MS_Model_Fluid3D::BuildSystem( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::BuildSystem() \n";
#endif

    //Build constant matrices
    M_Fluid->buildSystem();

    //Initialize problem coefficients
    M_alpha  = 0.0;
    *M_beta *= 0.0;
    *M_RHS  *= 0.0;

    //Initialize velocity and pressure to zero
    M_Fluid->initialize( *M_FluidSolution );

    //Set problem coefficients
    M_Fluid->updateSystem( M_alpha, *M_beta, *M_RHS );

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Model_Fluid3D::UpdateSystem( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::UpdateSystem() \n";
#endif

    //BDF
    if ( M_FluidData->isFirstTimeStep() )
        M_FluidBDF->bdf_u().initialize_unk( M_Fluid->solution() );
    else
        M_FluidBDF->bdf_u().shift_right( M_Fluid->solution() );

    //Time
    M_FluidData->updateTime();

    //Problem coefficients
    M_alpha = M_FluidBDF->bdf_u().coeff_der( 0 ) / M_FluidData->getTimeStep();
    *M_beta = M_FluidBDF->bdf_u().extrap();
    *M_RHS  = M_Fluid->matrMass() * M_FluidBDF->bdf_u().time_der( M_FluidData->getTimeStep() );

    //Set problem coefficients
    M_Fluid->updateSystem( M_alpha, *M_beta, *M_RHS );

    //Recompute preconditioner
    M_Fluid->resetPrec( true ); //TODO: Recompute only each N time step!

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Model_Fluid3D::SolveSystem( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SolveSystem() \n";
#endif

    //Solve the problem
    M_Fluid->iterate( *M_FluidBC->Handler_ptr() );

    UInt subITMax  = 0;   //TODO: Move this on a file!
    Real tolerance = 1.e-6; //TODO: Move this on a file!
    Real residual;
    for ( UInt subIT = 0; subIT < subITMax; ++subIT )
    {
        residual = ( *M_beta - M_Fluid->solution() ).Norm2(); // Residual is computed on the whole solution vector

        if ( M_displayer->isLeader() )
        {
            std::cout << "  F-  Sub-iteration n.:                        " << subIT << std::endl;
            std::cout << "  F-  Residual:                                " << residual << std::endl;
        }

        // Verify tolerance
        if ( residual <= tolerance )
            break;

        //Picard iteration for NonLinear Navier-Stokes
        *M_beta = M_Fluid->solution();
        M_Fluid->updateSystem( M_alpha, *M_beta, *M_RHS );

        M_Fluid->iterate( *M_FluidBC->Handler_ptr() );
    }

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Model_Fluid3D::SaveSolution( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SaveSolution() \n";
#endif

    //Post-processing
    *M_FluidSolution = M_Fluid->solution();
    M_output->postProcess( M_FluidData->getTime() );

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Model_Fluid3D::ShowMe( void )
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "uOrder             = " << M_FluidData->uOrder() << std::endl
                  << "pOrder             = " << M_FluidData->pOrder() << std::endl << std::endl;

        std::cout << "uDOF               = " << 3 * M_uDOF << std::endl
                  << "pDOF               = " << M_pDOF << std::endl << std::endl << std::endl << std::endl;
    }

    //MPI Barrier
    M_comm->Barrier();
}

// ===================================================
//! Private Methods
// ===================================================
void
MS_Model_Fluid3D::setupFEspace( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::setupFEspace() \n";
#endif

    //Transform mesh
    M_FluidData->mesh()->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    //Partition mesh
    M_FluidMesh.reset( new FluidMeshType( *M_FluidData->mesh(), *M_comm ) );
    M_FluidData->setMesh( M_FluidMesh->mesh() );

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

    M_uFESpace.reset( new FESpaceType( *M_FluidMesh, *u_refFE, *u_qR, *u_bdQr, 3, *M_comm ) );
    M_pFESpace.reset( new FESpaceType( *M_FluidMesh, *p_refFE, *p_qR, *p_bdQr, 1, *M_comm ) );
}

void
MS_Model_Fluid3D::setupDOF( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::setupDOF \n";
#endif

    //DOF
    M_uDOF = M_uFESpace->dof().numTotalDof();
    M_pDOF = M_pFESpace->dof().numTotalDof();
    //M_uDOF = M_uFESpace->map().getMap(Unique)->NumGlobalElements();
    //M_pDOF = M_pFESpace->map().getMap(Unique)->NumGlobalElements();
}

UInt
MS_Model_Fluid3D::imposeFluxes( const boost::shared_ptr< FluidBCType >& BC )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::imposeFluxes \n";
#endif

    std::vector< BCName > FluxVector = BC->Handler_ptr()->getBCWithType( Flux );

    UInt offset = M_uFESpace->map().getMap( Unique )->NumGlobalElements() + M_pFESpace->map().getMap( Unique )->NumGlobalElements();
    UInt numLM = static_cast< UInt > ( FluxVector.size() );
    for ( UInt i = 0; i < numLM; ++i )
        BC->Handler_ptr()->setOffset( FluxVector[i], offset + i );

    return numLM;
}

} // Namespace LifeV
