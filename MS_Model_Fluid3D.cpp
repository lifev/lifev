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
    M_fluid         (),
    M_fluidBC       (),
    M_fluidBDF      (),
    M_fluidData     (),
    M_fluidMesh     (),
    M_fluidFullMap  (),
    M_fluidSolution (),
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
    M_fluid         ( Fluid3D.M_fluid ),
    M_fluidBC       ( Fluid3D.M_fluidBC ),
    M_fluidBDF      ( Fluid3D.M_fluidBDF ),
    M_fluidData     ( Fluid3D.M_fluidData ),
    M_fluidMesh     ( Fluid3D.M_fluidMesh ),
    M_fluidFullMap  ( Fluid3D.M_fluidFullMap ),
    M_fluidSolution ( Fluid3D.M_fluidSolution ),
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
MS_Model_Fluid3D::operator=( const MS_Model_Fluid3D& fluid3D )
{
    if ( this != &fluid3D )
    {
        super::operator=( fluid3D );
        M_output        = fluid3D.M_output;
        M_fluid         = fluid3D.M_fluid;
        M_fluidBC       = fluid3D.M_fluidBC;
        M_fluidBDF      = fluid3D.M_fluidBDF;
        M_fluidMesh     = fluid3D.M_fluidMesh;
        M_fluidFullMap  = fluid3D.M_fluidFullMap;
        M_fluidSolution = fluid3D.M_fluidSolution;
        M_uFESpace      = fluid3D.M_uFESpace;
        M_pFESpace      = fluid3D.M_pFESpace;
        M_uDOF          = fluid3D.M_uDOF;
        M_pDOF          = fluid3D.M_pDOF;
        M_alpha         = fluid3D.M_alpha;
        M_beta          = fluid3D.M_beta;
        M_RHS           = fluid3D.M_RHS;
    }

    return *this;
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
    M_fluidData.reset( new fluidData_type( M_dataFile ) );

    //dataPhysics
    if ( !M_dataFile.checkVariable( "fluid/physics/density" ) )
        M_fluidData->density( M_dataPhysics->GetFluidDensity() );
    if ( !M_dataFile.checkVariable( "fluid/physics/viscosity" ) )
        M_fluidData->viscosity( M_dataPhysics->GetFluidViscosity() );

    //dataTime
    M_fluidData->setInitialTime( M_dataTime->getInitialTime() );
    M_fluidData->setEndTime( M_dataTime->getEndTime() );

    M_fluidData->setTime( M_dataTime->getInitialTime() );
    M_fluidData->setTimeStep( M_dataTime->getTimeStep() );
    //M_fluidData->setTimeStep	( globalTimeStep / std::ceil( globalTimeStep / M_fluidData->getTimeStep() ) );

    //FEspace & DOF
    setupFEspace();
    setupDOF();

    //Boundary Conditions
    M_fluidBC.reset( new fluidBC_type( M_dataFile, "fluid" ) );
    M_fluidBC->SetOperator( M_fluid ); //MUST BE MOVED AFTER M_fluid.reset !!!
    M_fluidBC->BuildHandler();
    M_fluidBC->UpdateOperatorVariables(); //MUST BE MOVED INSIDE THE UPDATE !!!
}

void
MS_Model_Fluid3D::SetupModel( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupProblem() \n";
#endif

    //Fluid
    M_fluid.reset( new fluid_type( *M_fluidData, *M_uFESpace, *M_pFESpace, *M_comm, imposeFluxes() ) );

    M_fluid->setUp( M_dataFile ); //Remove Preconditioner and Solver if possible!

    //Fluid MAP
    M_fluidFullMap.reset( new EpetraMap( M_fluid->getMap() ) );

    //BDF
    M_fluidBDF.reset( new fluidBDF_type( M_fluidData->getBDF_order() ) );

    //Problem coefficients
    M_beta.reset( new fluidVector_type( M_fluidFullMap ) );
    M_RHS.reset ( new fluidVector_type( M_fluidFullMap ) );

    //Post-processing
    M_output.reset( new output_type( M_dataFile, M_fluidMesh->mesh(), M_modelName + "_" + number2string( M_ID ), M_comm->MyPID() ) );

    M_fluidSolution.reset( new fluidVector_type( M_fluidFullMap, Repeated ) );
    //M_fluidSolution.reset( new fluidVector_type( M_fluid->solution(), Repeated ) );

    M_output->addVariable( ExporterData::Vector, "velocity", M_fluidSolution, UInt( 0 ), M_uDOF );
    M_output->addVariable( ExporterData::Scalar, "pressure", M_fluidSolution, 3 * M_uDOF, 3 * M_uDOF + M_pDOF );

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Model_Fluid3D::BuildSystem( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::BuildSystem() \n";
#endif

    //Build constant matrices
    M_fluid->buildSystem();

    //Initialize problem coefficients
    M_alpha  = 0.0;
    *M_beta *= 0.0;
    *M_RHS  *= 0.0;

    //Initialize velocity and pressure to zero
    M_fluid->initialize( *M_fluidSolution );

    //Set problem coefficients
    M_fluid->updateSystem( M_alpha, *M_beta, *M_RHS );

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
    if ( M_fluidData->isFirstTimeStep() )
        M_fluidBDF->bdf_u().initialize_unk( M_fluid->solution() );
    else
        M_fluidBDF->bdf_u().shift_right( M_fluid->solution() );

    //Time
    M_fluidData->updateTime();

    //Problem coefficients
    M_alpha = M_fluidBDF->bdf_u().coeff_der( 0 ) / M_fluidData->getTimeStep();
    *M_beta = M_fluidBDF->bdf_u().extrap();
    *M_RHS  = M_fluid->matrMass() * M_fluidBDF->bdf_u().time_der( M_fluidData->getTimeStep() );

    //Set problem coefficients
    M_fluid->updateSystem( M_alpha, *M_beta, *M_RHS );

    //Recompute preconditioner
    M_fluid->resetPrec( true ); // TODO: Recompute only each N time step!

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
    M_fluid->iterate( *M_fluidBC->Handler_ptr() );

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
    *M_fluidSolution = M_fluid->solution();
    M_output->postProcess( M_fluidData->getTime() );

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Model_Fluid3D::ShowMe( void )
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "uOrder             = " << M_fluidData->uOrder() << std::endl
                  << "pOrder             = " << M_fluidData->pOrder() << std::endl << std::endl;

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
    M_fluidData->mesh()->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    //Partition mesh
    M_fluidMesh.reset( new fluidMesh_type( *M_fluidData->mesh(), *M_comm ) );
    M_fluidData->setMesh( M_fluidMesh->mesh() );

    //Velocity FE Space
    const RefFE* u_refFE;
    const QuadRule* u_qR;
    const QuadRule* u_bdQr;

    if ( M_fluidData->uOrder().compare( "P2" ) == 0 )
    {
        u_refFE = &feTetraP2;
        u_qR = &quadRuleTetra15pt; // DoE 5
        u_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else
        if ( M_fluidData->uOrder().compare( "P1" ) == 0 )
        {
            u_refFE = &feTetraP1;
            u_qR = &quadRuleTetra4pt; // DoE 2
            u_bdQr = &quadRuleTria3pt; // DoE 2
        }
        else
            if ( M_fluidData->uOrder().compare( "P1Bubble" ) == 0 )
            {
                u_refFE = &feTetraP1bubble;
                u_qR = &quadRuleTetra64pt; // DoE 2
                u_bdQr = &quadRuleTria3pt; // DoE 2
            }
            else
            {
                if ( M_displayer->isLeader() )
                    std::cout << M_fluidData->uOrder() << " Velocity FE not implemented yet." << std::endl;
                exit( EXIT_FAILURE );
            }

    //Pressure FE Space
    const RefFE* p_refFE;
    const QuadRule* p_qR;
    const QuadRule* p_bdQr;

    if ( M_fluidData->pOrder().compare( "P2" ) == 0 )
    {
        p_refFE = &feTetraP2;
        p_qR = u_qR;
        p_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else
        if ( M_fluidData->pOrder().compare( "P1" ) == 0 )
        {
            p_refFE = &feTetraP1;
            p_qR = u_qR;
            p_bdQr = &quadRuleTria3pt; // DoE 2
        }
        else
        {
            if ( M_displayer->isLeader() )
                std::cout << M_fluidData->pOrder() << " pressure FE not implemented yet." << std::endl;
            exit( EXIT_FAILURE );
        }

    M_uFESpace.reset( new FESpace_type( *M_fluidMesh, *u_refFE, *u_qR, *u_bdQr, 3, *M_comm ) );
    M_pFESpace.reset( new FESpace_type( *M_fluidMesh, *p_refFE, *p_qR, *p_bdQr, 1, *M_comm ) );
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
MS_Model_Fluid3D::imposeFluxes( void )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::imposeFluxes \n";
#endif

    std::vector< BCName > fluxVector = M_fluidBC->Handler_ptr()->getBCWithType( Flux );
    UInt numLM = static_cast< UInt > ( fluxVector.size() );

    UInt offset = M_uFESpace->map().getMap( Unique )->NumGlobalElements() + M_pFESpace->map().getMap( Unique )->NumGlobalElements();

    for ( UInt i = 0; i < numLM; ++i )
        M_fluidBC->Handler_ptr()->setOffset( fluxVector[i], offset + i );

    return numLM;
}

} // Namespace LifeV
