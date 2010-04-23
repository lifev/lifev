//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

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
 *  @brief MultiScale Model 1D
 *
 *  @version 1.0
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *  @date 02-25-2010
 *
 *  @version 1.1
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 14-04-2010
 */

#include "MS_Model_1D.hpp"

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Model_1D::MS_Model_1D() :
    super                          (),
    M_Output                       (),
    M_Data                         ( new Data_Type() ),
    M_Physics                      (),
    M_Flux                         (),
    M_Source                       (),
    M_Solver                       (),
    M_BC                           (),
    M_Solution                     (),
    M_FESpace                      (),
    M_LinearSolver                 ()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_1D::MS_Model_1D() \n";
#endif

    //M_type = oneD;

    //Define the maps of the OneDimensionalModel objects
    OneDimensionalModel_MapsDefinition();

    //Register the objects
    Factory_OneDimensionalModel_Physics::instance().registerProduct( OneD_LinearPhysics,    &Create_OneDimensionalModel_Physics_Linear );
    Factory_OneDimensionalModel_Physics::instance().registerProduct( OneD_NonLinearPhysics, &Create_OneDimensionalModel_Physics_NonLinear );

    Factory_OneDimensionalModel_Flux::instance().registerProduct(    OneD_LinearFlux,       &Create_OneDimensionalModel_Flux_Linear );
    Factory_OneDimensionalModel_Flux::instance().registerProduct(    OneD_NonLinearFlux,    &Create_OneDimensionalModel_Flux_NonLinear );

    Factory_OneDimensionalModel_Source::instance().registerProduct(  OneD_LinearSource,     &Create_OneDimensionalModel_Source_Linear );
    Factory_OneDimensionalModel_Source::instance().registerProduct(  OneD_NonLinearSource,  &Create_OneDimensionalModel_Source_NonLinear );
}

// ===================================================
// MultiScale PhysicalModel Virtual Methods
// ===================================================
void
MS_Model_1D::SetupData( const std::string& FileName )
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupData( ) \n";
#endif

    super::SetupData( FileName );

    GetPot DataFile( FileName );

    M_Data->setup( DataFile );
    //M_Data->showMe();

    M_LinearSolver.reset( new LinearSolver_Type( *M_comm ) );
    M_LinearSolver->setUpPrec        ( DataFile, "1D_Model/prec" );
    M_LinearSolver->setDataFromGetPot( DataFile, "1D_Model/solver");
}

void
MS_Model_1D::SetupModel()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupProblem() \n";
#endif

    M_Physics = Physics_PtrType( Factory_OneDimensionalModel_Physics::instance().createObject( M_Data->PhysicsType() ) );
    M_Physics->SetData( M_Data );

    M_Flux = Flux_PtrType( Factory_OneDimensionalModel_Flux::instance().createObject( M_Data->FluxType() ) );
    M_Flux->SetPhysics( M_Physics );

    M_Source = Source_PtrType( Factory_OneDimensionalModel_Source::instance().createObject( M_Data->SourceType() ) );
    M_Source->SetPhysics( M_Physics );

    SetupFESpace();

    M_Solver.reset ( new Solver_Type ( M_Physics, M_Flux, M_Source, M_FESpace, M_comm ) );
    M_Solver->setup();

    M_Solver->setLinearSolver( M_LinearSolver );

    M_BC.reset( new BC_Type( M_Solver->U_thistime(), M_Flux, M_FESpace->dim() ) );
}

void
MS_Model_1D::BuildSystem()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_1D::BuildSystem() \n";
#endif

    M_Solver->initialize();
}

void
MS_Model_1D::UpdateSystem()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_1D::UpdateSystem() \n";
#endif

    M_Solver->timeAdvance();
}

void
MS_Model_1D::SolveSystem()
{
#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_1D::SolveSystem() \n";
#endif

    M_Solver->iterate( *M_BC, M_Data->dataTime()->getTime() );
}

void
MS_Model_1D::SaveSolution()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_1D::SaveSolution() \n";
#endif

    //Post-processing

//     *M_FluidSolution = M_Fluid->solution();
//     M_Output->postProcess( M_FluidData->getTime() );

}

void
MS_Model_1D::ShowMe()
{
}

// ===================================================
// Methods
// ===================================================
void
MS_Model_1D::SetupLinearData( const std::string& /*FileName*/ )
{

}

void
MS_Model_1D::SetupLinearModel()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupLinearModel( ) \n";
#endif
}

void
MS_Model_1D::UpdateLinearModel()
{
}

void
MS_Model_1D::SolveLinearModel( bool& /*SolveLinearSystem*/ )
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_1D::SolveLinearModel() \n";
#endif

}

// ===================================================
// Get Methods
// ===================================================
MS_Model_1D::BC_Type&
MS_Model_1D::GetBC()
{
    return *M_BC;
}

MS_Model_1D::Data_Type&
MS_Model_1D::GetData()
{
    return *M_Data;
}

MS_Model_1D::Physics_PtrType&
MS_Model_1D::GetPhysics()
{
    return M_Physics;
}

MS_Model_1D::Flux_PtrType&
MS_Model_1D::GetFlux()
{
    return M_Flux;
}

MS_Model_1D::Source_PtrType&
MS_Model_1D::GetSource()
{
    return M_Source;
}

MS_Model_1D::FESpace_Type&
MS_Model_1D::GetFESpace()
{
    return *M_FESpace;
}

MS_Model_1D::Solver_Type&
MS_Model_1D::GetSolver()
{
    return *M_Solver;
}

// ===================================================
// Set Methods
// ===================================================
void
MS_Model_1D::SetBC( boost::shared_ptr<BC_Type>& BC )
{
    M_BC = BC;
}

// ===================================================
// Private Methods
// ===================================================
void
MS_Model_1D::SetupFESpace()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_1D::setupFEspace() \n";
#endif

    const RefFE*    refFE = &feSegP1;
    const QuadRule* qR    = &quadRuleSeg3pt;
    const QuadRule* bdQr  = &quadRuleSeg1pt;

    M_FESpace.reset( new FESpace_Type( M_Data->mesh(), *refFE, *qR, *bdQr, 1, *M_comm ) );
}


} // Namespace LifeV
