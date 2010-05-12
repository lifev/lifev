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

#include "MS_Model_FSI1D.hpp"

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Model_FSI1D::MS_Model_FSI1D() :
    super                          (),
#ifdef HAVE_HDF5
    M_Exporter                     ( new IOFile_Type() ),
    M_ExporterMesh                 ( new Mesh_Type() ),
#endif
    M_Data                         ( new Data_Type() ),
    M_Physics                      (),
    M_Flux                         (),
    M_Source                       (),
    M_Solver                       ( new Solver_Type() ),
    M_BC                           ( new BCInterface_Type() ),
    M_Solution                     (),
    M_FESpace                      (),
    M_LinearSolver                 ()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_FSI1D::MS_Model_FSI1D() \n";
#endif

    M_type = FSI1D;

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

MS_Model_FSI1D::MS_Model_FSI1D( const MS_Model_FSI1D& FSI1D ) :
    super                          ( FSI1D ),
#ifdef HAVE_HDF5
    M_Exporter                     ( FSI1D.M_Exporter ),
    M_ExporterMesh                 ( FSI1D.M_ExporterMesh ),
#endif
    M_Data                         ( FSI1D.M_Data ),
    M_Physics                      ( FSI1D.M_Physics ),
    M_Flux                         ( FSI1D.M_Flux ),
    M_Source                       ( FSI1D.M_Source ),
    M_Solver                       ( FSI1D.M_Solver ),
    M_BC                           ( FSI1D.M_BC ),
    M_Solution                     ( FSI1D.M_Solution ),
    M_FESpace                      ( FSI1D.M_FESpace ),
    M_LinearSolver                 ( FSI1D.M_LinearSolver )
{
}

// ===================================================
// Operators
// ===================================================
MS_Model_FSI1D&
MS_Model_FSI1D::operator=( const MS_Model_FSI1D& FSI1D )
{
    if ( this != &FSI1D )
    {
        super::operator=( FSI1D );
#ifdef HAVE_HDF5
        M_Exporter                     = FSI1D.M_Exporter;
        M_ExporterMesh                 = FSI1D.M_ExporterMesh;
#endif
        M_Data                         = FSI1D.M_Data;
        M_Physics                      = FSI1D.M_Physics;
        M_Flux                         = FSI1D.M_Flux;
        M_Source                       = FSI1D.M_Source;
        M_Solver                       = FSI1D.M_Solver;
        M_BC                           = FSI1D.M_BC;
        M_Solution                     = FSI1D.M_Solution;
        M_FESpace                      = FSI1D.M_FESpace;
        M_LinearSolver                 = FSI1D.M_LinearSolver;
    }

    return *this;
}

// ===================================================
// MultiScale PhysicalModel Virtual Methods
// ===================================================
void
MS_Model_FSI1D::SetupData( const std::string& FileName )
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_FSI1D::SetupData( ) \n";
#endif

    super::SetupData( FileName );

    GetPot DataFile( FileName );

    M_Data->setup( DataFile );

    //If global physical quantities are specified also in the data file, replace them
    if ( !DataFile.checkVariable( "1D_Model/PhysicalParameters/density" ) )
        M_Data->setDensity( M_dataPhysics->GetFluidDensity() );
    if ( !DataFile.checkVariable( "1D_Model/PhysicalParameters/viscosity" ) )
        M_Data->setViscosity( M_dataPhysics->GetFluidViscosity() );
    if ( !DataFile.checkVariable( "1D_Model/PhysicalParameters/poisson" ) )
        M_Data->setPoisson( M_dataPhysics->GetStructurePoissonCoefficient() );
    if ( !DataFile.checkVariable( "1D_Model/PhysicalParameters/young" ) )
        M_Data->setYoung( M_dataPhysics->GetStructureYoungModulus() );

    //After changing some parameters we need to update the coefficients
    M_Data->UpdateCoefficients();

    //1D Model Physics
    M_Physics = Physics_PtrType( Factory_OneDimensionalModel_Physics::instance().createObject( M_Data->PhysicsType() ) );
    M_Physics->SetData( M_Data );

    //1D Model Flux
    M_Flux = Flux_PtrType( Factory_OneDimensionalModel_Flux::instance().createObject( M_Data->FluxType() ) );
    M_Flux->SetPhysics( M_Physics );

    //1D Model Source
    M_Source = Source_PtrType( Factory_OneDimensionalModel_Source::instance().createObject( M_Data->SourceType() ) );
    M_Source->SetPhysics( M_Physics );

    //Linear Solver
    M_LinearSolver.reset( new LinearSolver_Type( *M_comm ) );
    M_LinearSolver->setUpPrec        ( DataFile, "1D_Model/prec" );
    M_LinearSolver->setDataFromGetPot( DataFile, "1D_Model/solver");

    //1D Model Solver
    M_Solver->setCommunicator( M_comm );
    M_Solver->setProblem( M_Physics, M_Flux, M_Source );
    M_Solver->setLinearSolver( M_LinearSolver );

    //BC - We need to create the BCHandler before using it
    M_BC->SetOperator( M_Solver );
    M_BC->CreateHandler();
    //M_BC->FillHandler( FileName, "1D_Model" );

#ifdef HAVE_HDF5
    //Exporter
    M_Exporter->setDataFromGetPot( DataFile );
    M_Exporter->setPrefix( "Step_" + number2string( MS_ProblemStep ) + "_Model_" + number2string( M_ID ) );
    M_Exporter->setDirectory( MS_ProblemFolder );

    M_ExporterMesh->setup( M_Data->Length(), M_Data->nbElem() );
#endif

}

void
MS_Model_FSI1D::SetupModel()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_FSI1D::SetupProblem() \n";
#endif

    //FEspace
    SetupFESpace();

    //1D Model Solver setup
    M_Solver->setup();

    //Set default BC (has to be called after setting other BC)
    M_BC->GetHandler()->setDefaultBC( M_Flux, M_Source, M_Solver->U_thistime() );

#ifdef HAVE_HDF5
    //Post-processing
    M_Exporter->setMeshProcId( M_ExporterMesh, M_comm->MyPID() );
    M_Solution.resize( 5 );
    for ( UInt i(0) ; i < static_cast <UInt> ( M_Solution.size() ) ; ++i)
    {
        M_Solution[i].reset( new Vector_Type( (*M_Solver->U_thistime())[i], M_Exporter->mapType() ) );
        if ( M_Exporter->mapType() == Unique )
            M_Solution[i]->setCombineMode( Zero );
    }

    M_Exporter->addVariable( ExporterData::Scalar, "Area",      M_Solution[0], static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_Exporter->addVariable( ExporterData::Scalar, "Flow Rate", M_Solution[1], static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_Exporter->addVariable( ExporterData::Scalar, "W1",        M_Solution[2], static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_Exporter->addVariable( ExporterData::Scalar, "W2",        M_Solution[3], static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_Exporter->addVariable( ExporterData::Scalar, "Pressure",  M_Solution[4], static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
#endif
}

void
MS_Model_FSI1D::BuildSystem()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_FSI1D::BuildSystem() \n";
#endif

    M_Solver->initialize();
}

void
MS_Model_FSI1D::UpdateSystem()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_FSI1D::UpdateSystem() \n";
#endif

    M_Solver->timeAdvance();
}

void
MS_Model_FSI1D::SolveSystem()
{
#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_FSI1D::SolveSystem() \n";
#endif

    M_BC->UpdateOperatorVariables();
    M_Solver->iterate( *M_BC->GetHandler() );
}

void
MS_Model_FSI1D::SaveSolution()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_FSI1D::SaveSolution() \n";
#endif

#ifdef HAVE_HDF5
    //Post-processing
    for ( UInt i(0) ; i < static_cast <UInt> ( M_Solution.size() ) ; ++i)
        *(M_Solution[i]) = (*M_Solver->U_thistime())[i];

    M_Exporter->postProcess( M_Data->dataTime()->getTime() );

    if ( M_Data->dataTime()->isLastTimeStep() )
        M_Exporter->CloseFile();
#endif

}

void
MS_Model_FSI1D::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "FE space            = " << "P1" << std::endl
                  << "DOF                 = " << M_Data->mesh()->numPoints() << std::endl << std::endl;

        std::cout << "maxH                = " << M_Data->mesh()->maxH() << std::endl
                  << "meanH               = " << M_Data->mesh()->meanH() << std::endl << std::endl;
    }
}

// ===================================================
// Methods
// ===================================================
void
MS_Model_FSI1D::SetupLinearData( const std::string& /*FileName*/ )
{

}

void
MS_Model_FSI1D::SetupLinearModel()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_FSI1D::SetupLinearModel( ) \n";
#endif
}

void
MS_Model_FSI1D::UpdateLinearModel()
{
}

void
MS_Model_FSI1D::SolveLinearModel( bool& /*SolveLinearSystem*/ )
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_FSI1D::SolveLinearModel() \n";
#endif

}

// ===================================================
// Get Methods (couplings)
// ===================================================
MS_Model_FSI1D::BCInterface_Type&
MS_Model_FSI1D::GetBCInterface() const
{
    return *M_BC;
}

Real
MS_Model_FSI1D::GetBoundaryDensity( const BCFlag& /*Flag*/ ) const
{
    return M_Data->DensityRho();
}

Real
MS_Model_FSI1D::GetBoundaryViscosity( const BCFlag& /*Flag*/ ) const
{
    return M_Data->Viscosity();
}

Real
MS_Model_FSI1D::GetBoundaryArea( const BCFlag& Flag ) const
{
    return M_Solver->BoundaryValue( OneD_A, (Flag == 0) ? OneD_left : OneD_right );
}

Real
MS_Model_FSI1D::GetBoundaryFlowRate( const BCFlag& Flag ) const
{
    return M_Solver->BoundaryValue( OneD_Q, (Flag == 0) ? OneD_left : OneD_right );
}

Real
MS_Model_FSI1D::GetBoundaryPressure( const BCFlag& Flag ) const
{
    return M_Solver->BoundaryValue( OneD_P, (Flag == 0) ? OneD_left : OneD_right );
}

Real
MS_Model_FSI1D::GetBoundaryDynamicPressure( const BCFlag& Flag ) const
{
    return 0.5 * GetBoundaryDensity( Flag ) * ( GetBoundaryFlowRate( Flag ) * GetBoundaryFlowRate( Flag ) )
                                            / ( GetBoundaryArea( Flag ) * GetBoundaryArea( Flag ) );
}

Real
MS_Model_FSI1D::GetBoundaryStress( const BCFlag& Flag, const stressTypes& StressType ) const
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

        default:

            std::cout << "ERROR: Invalid stress type [" << Enum2String( StressType, stressMap ) << "]" << std::endl;

            return 0.0;
    }
}

// ===================================================
// Get Methods
// ===================================================
MS_Model_FSI1D::BC_Type&
MS_Model_FSI1D::GetBC() const
{
    return *(M_BC->GetHandler());
}


MS_Model_FSI1D::Data_Type&
MS_Model_FSI1D::GetData() const
{
    return *M_Data;
}

MS_Model_FSI1D::Physics_PtrType
MS_Model_FSI1D::GetPhysics() const
{
    return M_Physics;
}

MS_Model_FSI1D::Flux_PtrType
MS_Model_FSI1D::GetFlux() const
{
    return M_Flux;
}

MS_Model_FSI1D::Source_PtrType
MS_Model_FSI1D::GetSource() const
{
    return M_Source;
}

MS_Model_FSI1D::FESpace_PtrType
MS_Model_FSI1D::GetFESpace() const
{
    return M_FESpace;
}

MS_Model_FSI1D::Solver_PtrType
MS_Model_FSI1D::GetSolver() const
{
    return M_Solver;
}

const MS_Model_FSI1D::Solution_PtrType
MS_Model_FSI1D::GetSolution() const
{
    return M_Solver->U_thistime();
}

// ===================================================
// Private Methods
// ===================================================
void
MS_Model_FSI1D::SetupFESpace()
{

#ifdef DEBUG
    Debug( 8130 ) << "MS_Model_FSI1D::setupFEspace() \n";
#endif

    //Transform mesh
    boost::array< Real, NDIM > NullTransformation;
    NullTransformation[0] = 0.; NullTransformation[1] = 0.; NullTransformation[2] = 0.;

    //The real mesh can be only scaled due to OneDimensionalModel_Solver conventions
    M_Data->mesh()->transformMesh( M_geometryScale, NullTransformation, NullTransformation );

#ifdef HAVE_HDF5
    //The mesh for the post-processing can be rotated
    M_ExporterMesh->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );
#endif

    //Setup FESpace
    const RefFE*    refFE = &feSegP1;
    const QuadRule* qR    = &quadRuleSeg3pt;
    const QuadRule* bdQr  = &quadRuleSeg1pt;

    M_FESpace.reset( new FESpace_Type( M_Data->mesh(), *refFE, *qR, *bdQr, 1, *M_comm ) );
    M_Solver->setFESpace( M_FESpace );
}

} // Namespace LifeV
