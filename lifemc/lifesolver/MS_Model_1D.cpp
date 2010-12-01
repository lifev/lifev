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
 *  @version 1.1
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *  @date 26-02-2010
 *
 *  @version 1.2 and subsequents
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 23-04-2010
 */

#include "MS_Model_1D.hpp"

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Model_1D::MS_Model_1D() :
        super                          (),
#ifdef HAVE_HDF5
        M_Exporter                     ( new IOFile_Type() ),
        M_Importer                     ( new IOFile_Type() ),
        M_ExporterMesh                 ( new Mesh_Type() ),
#endif
#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
        M_LinearBC                     ( new BC_Type() ),
        M_LinearSolution               ( new Solution_Type() ),
        M_BCPreviousTimeSteps          (),
        M_BCBaseDelta                  (),
        M_BCDelta                      (),
        M_BCDeltaType                  (),
        M_BCDeltaSide                  (),
#endif
        M_Data                         ( new Data_Type() ),
        M_BC                           ( new BCInterface_Type() ),
        M_Physics                      (),
        M_Flux                         (),
        M_Source                       (),
        M_Solver                       ( new Solver_Type() ),
        M_LinearSolver                 (),
        M_FESpace                      (),
        M_Solution_tn                  ( new Solution_Type() ),
        M_Solution                     ( new Solution_Type() )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::MS_Model_1D() \n";
#endif

    M_type = OneDimensional;

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

MS_Model_1D::MS_Model_1D( const MS_Model_1D& OneDimensionalModel ) :
        super                          ( OneDimensionalModel ),
#ifdef HAVE_HDF5
        M_Exporter                     ( OneDimensionalModel.M_Exporter ),
        M_Importer                     ( OneDimensionalModel.M_Importer ),
        M_ExporterMesh                 ( OneDimensionalModel.M_ExporterMesh ),
#endif
#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
        M_LinearBC                     ( OneDimensionalModel.M_LinearBC ),
        M_LinearSolution               ( OneDimensionalModel.M_LinearSolution ),
        M_BCPreviousTimeSteps          ( OneDimensionalModel.M_BCPreviousTimeSteps ),
        M_BCBaseDelta                  ( OneDimensionalModel.M_BCBaseDelta ),
        M_BCDelta                      ( OneDimensionalModel.M_BCDelta ),
        M_BCDeltaType                  ( OneDimensionalModel.M_BCDeltaType ),
        M_BCDeltaSide                  ( OneDimensionalModel.M_BCDeltaSide ),
#endif
        M_Data                         ( OneDimensionalModel.M_Data ),
        M_BC                           ( OneDimensionalModel.M_BC ),
        M_Physics                      ( OneDimensionalModel.M_Physics ),
        M_Flux                         ( OneDimensionalModel.M_Flux ),
        M_Source                       ( OneDimensionalModel.M_Source ),
        M_Solver                       ( OneDimensionalModel.M_Solver ),
        M_LinearSolver                 ( OneDimensionalModel.M_LinearSolver ),
        M_FESpace                      ( OneDimensionalModel.M_FESpace ),
        M_Solution_tn                  ( OneDimensionalModel.M_Solution_tn ),
        M_Solution                     ( OneDimensionalModel.M_Solution )
{
}

// ===================================================
// Operators
// ===================================================
MS_Model_1D&
MS_Model_1D::operator=( const MS_Model_1D& OneDimensionalModel )
{
    if ( this != &OneDimensionalModel )
    {
        super::operator=( OneDimensionalModel );
#ifdef HAVE_HDF5
        M_Exporter                     = OneDimensionalModel.M_Exporter;
        M_Importer                     = OneDimensionalModel.M_Importer;
        M_ExporterMesh                 = OneDimensionalModel.M_ExporterMesh;
#endif
#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
        M_LinearBC                     = OneDimensionalModel.M_LinearBC;
        M_LinearSolution               = OneDimensionalModel.M_LinearSolution;
        M_BCPreviousTimeSteps          = OneDimensionalModel.M_BCPreviousTimeSteps;
        M_BCBaseDelta                  = OneDimensionalModel.M_BCBaseDelta;
        M_BCDelta                      = OneDimensionalModel.M_BCDelta;
        M_BCDeltaType                  = OneDimensionalModel.M_BCDeltaType;
        M_BCDeltaSide                  = OneDimensionalModel.M_BCDeltaSide;
#endif
        M_Data                         = OneDimensionalModel.M_Data;
        M_BC                           = OneDimensionalModel.M_BC;
        M_Physics                      = OneDimensionalModel.M_Physics;
        M_Flux                         = OneDimensionalModel.M_Flux;
        M_Source                       = OneDimensionalModel.M_Source;
        M_Solver                       = OneDimensionalModel.M_Solver;
        M_LinearSolver                 = OneDimensionalModel.M_LinearSolver;
        M_FESpace                      = OneDimensionalModel.M_FESpace;
        M_Solution_tn                  = OneDimensionalModel.M_Solution_tn;
        M_Solution                     = OneDimensionalModel.M_Solution;
    }

    return *this;
}

// ===================================================
// MultiScale PhysicalModel Virtual Methods
// ===================================================
void
MS_Model_1D::SetupData( const std::string& FileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupData( ) \n";
#endif

    // Preliminary setup of the communicator
#ifdef EPETRA_MPI
    MPI_Comm LocalComm;
    MPI_Comm_split( ( dynamic_cast<Epetra_MpiComm*> ( &(*M_comm) ) )->Comm(), M_comm->MyPID(), M_comm->MyPID(), &LocalComm );
    M_comm.reset( new Epetra_MpiComm( LocalComm ) );
#else
    M_comm.reset( new Epetra_SerialComm() );
#endif

    super::SetupData( FileName );

    GetPot DataFile( FileName );

    M_Data->setup( DataFile );
    if ( M_globalData.get() )
        SetupGlobalData( FileName );

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
    M_LinearSolver.reset( new LinearSolver_Type( M_comm ) );
    M_LinearSolver->setUpPrec        ( DataFile, "1D_Model/prec" );
    M_LinearSolver->setDataFromGetPot( DataFile, "1D_Model/solver" );
    M_LinearSolver->setParameters();

    //1D Model Solver
    M_Solver->setCommunicator( M_comm );
    M_Solver->setProblem( M_Physics, M_Flux, M_Source );
    M_Solver->setLinearSolver( M_LinearSolver );

    //BC - We need to create the BCHandler before using it
    M_BC->CreateHandler();
    //M_BC->FillHandler( FileName, "1D_Model" );

    //Exporters
    M_Data->setPostprocessingDirectory( MS_ProblemFolder );
    M_Data->setPostprocessingFile( "Step_" + number2string( MS_ProblemStep ) + "_Model_" + number2string( M_ID ) );

#ifdef HAVE_HDF5
    M_Exporter->setDataFromGetPot( DataFile );
    M_Exporter->setPrefix( "Step_" + number2string( MS_ProblemStep ) + "_Model_" + number2string( M_ID ) );
    M_Exporter->setDirectory( MS_ProblemFolder );

    M_ExporterMesh->setup( M_Data->Length(), M_Data->NumberOfElements() );

    M_Importer->setDataFromGetPot( DataFile );
    M_Importer->setPrefix( "Step_" + number2string( MS_ProblemStep - 1 ) + "_Model_" + number2string( M_ID ) );
    M_Importer->setDirectory( MS_ProblemFolder );
#endif

}

void
MS_Model_1D::SetupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupProblem() \n";
#endif

    //FEspace
    SetupFESpace();

    //Setup solution
    M_Solver->setupSolution( *M_Solution );
    M_Solver->setupSolution( *M_Solution_tn );

    //Set default BC (has to be called after setting other BC)
    M_BC->GetHandler()->setDefaultBC();
    M_BC->SetOperator( M_Solver );
    M_BC->SetSolution( M_Solution );
    M_BC->SetFluxSource( M_Flux, M_Source );

#ifdef HAVE_HDF5
    //Post-processing
    M_Exporter->setMeshProcId( M_ExporterMesh, M_comm->MyPID() );

    //M_Exporter->addVariable( ExporterData::Scalar, "Solid Area",      (*M_Solution)["A"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_Exporter->addVariable( ExporterData::Scalar, "Area ratio",      (*M_Solution)["A/A0-1"], static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_Exporter->addVariable( ExporterData::Scalar, "Fluid Flow Rate", (*M_Solution)["Q"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    //M_Exporter->addVariable( ExporterData::Scalar, "W1",              (*M_Solution)["W1"],   static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    //M_Exporter->addVariable( ExporterData::Scalar, "W2",              (*M_Solution)["W2"],   static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_Exporter->addVariable( ExporterData::Scalar, "Fluid Pressure",  (*M_Solution)["P"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
#endif

#ifdef HAVE_MATLAB_POSTPROCESSING
    //Matlab post-processing
    M_Solver->resetOutput( *M_Solution );
#endif

    //Setup solution
    InitializeSolution();

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
    {
        CreateLinearBC();
        UpdateLinearBC( *M_Solution );
        SetupLinearModel();
    }
#endif

}

void
MS_Model_1D::BuildSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::BuildSystem() \n";
#endif

    //M_Data->showMe();
    M_Solver->buildConstantMatrices();

    // Update previous solution
    UpdateSolution( *M_Solution, *M_Solution_tn );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
        UpdateLinearModel();
#endif

}

void
MS_Model_1D::UpdateSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::UpdateSystem() \n";
#endif

    // Update previous solution
    UpdateSolution( *M_Solution, *M_Solution_tn );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
        UpdateLinearModel();
#endif

}

void
MS_Model_1D::SolveSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SolveSystem() \n";
#endif

    Solve( *M_BC->GetHandler(), *M_Solution );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    if ( M_couplings.size() > 0 )
        UpdateLinearBC( *M_Solution );
#endif

}

void
MS_Model_1D::SaveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SaveSolution() \n";
#endif

#ifdef HAVE_HDF5
    M_Exporter->postProcess( M_Data->dataTime()->getTime() );

    if ( M_Data->dataTime()->isLastTimeStep() )
        M_Exporter->CloseFile();
#endif

#ifdef HAVE_MATLAB_POSTPROCESSING
    //Matlab post-processing
    M_Solver->postProcess( *M_Solution );
#endif

}

void
MS_Model_1D::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "FE order            = " << "P1" << std::endl
                  << "DOF                 = " << M_Data->mesh()->numPoints() << std::endl << std::endl;

        std::cout << "maxH                = " << M_Data->mesh()->maxH() << std::endl
                  << "meanH               = " << M_Data->mesh()->meanH() << std::endl << std::endl;
    }
}

// ===================================================
// Methods
// ===================================================
#ifdef JACOBIAN_WITH_FINITEDIFFERENCE

void
MS_Model_1D::SetupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupLinearModel( ) \n";
#endif

    // Define BCFunction for linear problem
    M_BCBaseDelta.setFunction( boost::bind( &MS_Model_1D::BCFunctionDelta, this, _1 ) );

    // The linear BCHandler is a copy of the original BCHandler with the LinearSolution instead of the true solution
    //M_LinearBC.reset( new BC_Type( *M_BC->GetHandler() ) ); // COPY CONSTRUCTOR NOT WORKING

    //Set left and right BC + default BC
    M_LinearBC->setBC( OneD_left, OneD_first, M_BC->GetHandler()->BC( OneD_left )->type( OneD_first ),
                       M_BC->GetHandler()->BC( OneD_left )->BCFunction( OneD_first ) );

    M_LinearBC->setBC( OneD_right, OneD_first, M_BC->GetHandler()->BC( OneD_right )->type( OneD_first ),
                       M_BC->GetHandler()->BC( OneD_right )->BCFunction( OneD_first ) );

    M_LinearBC->setDefaultBC();

    // Solution for the linear problem (this does not change anything in the solver)
    M_Solver->setupSolution( *M_LinearSolution );
    M_LinearBC->setSolution( M_LinearSolution );
    M_LinearBC->setFluxSource( M_Flux, M_Source );
}

void
MS_Model_1D::UpdateLinearModel()
{
    // The couplings should use the same value for the time interpolation order
    UInt timeInterpolationOrder( std::max( M_couplings[0]->GetTimeInterpolationOrder(), M_couplings[1]->GetTimeInterpolationOrder() ) );

    UInt containerSize( M_BCPreviousTimeSteps.size() );

    // If we have not yet enough samples for interpolation, we add a new one
    if ( containerSize <= timeInterpolationOrder )
    {
        ++containerSize;
        M_BCPreviousTimeSteps.push_back( M_BCPreviousTimeSteps[0] );
    }

    // Updating database
    for ( UInt i(1) ; i < containerSize ; ++i )
        M_BCPreviousTimeSteps[containerSize-i] = M_BCPreviousTimeSteps[containerSize-i-1];
}

void
MS_Model_1D::SolveLinearModel( bool& SolveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SolveLinearModel() \n";
#endif

    if ( !SolveLinearSystem )
        return;

    ImposePerturbation();

    Solve( *M_LinearBC, *M_LinearSolution, "L1D-" );

    ResetPerturbation();

    //This flag avoid recomputation of the same system
    SolveLinearSystem = false;
}

#endif

// ===================================================
// Get Methods (couplings)
// ===================================================
MS_Model_1D::BCInterface_Type&
MS_Model_1D::GetBCInterface() const
{
    return *M_BC;
}

Real
MS_Model_1D::GetBoundaryDensity( const BCFlag& /*Flag*/ ) const
{
    return M_Data->DensityRho();
}

Real
MS_Model_1D::GetBoundaryViscosity( const BCFlag& /*Flag*/ ) const
{
    return M_Data->Viscosity();
}

Real
MS_Model_1D::GetBoundaryArea( const BCFlag& Flag ) const
{
    return M_Solver->BoundaryValue( *M_Solution, OneD_A, FlagConverter( Flag ) );
}

Real
MS_Model_1D::GetBoundaryFlowRate( const BCFlag& Flag ) const
{
    return M_Solver->BoundaryValue( *M_Solution, OneD_Q, FlagConverter( Flag ) );
}

Real
MS_Model_1D::GetBoundaryPressure( const BCFlag& Flag ) const
{
    return M_Solver->BoundaryValue( *M_Solution, OneD_P, FlagConverter( Flag ) );
}

Real
MS_Model_1D::GetBoundaryDynamicPressure( const BCFlag& Flag ) const
{
    return 0.5 * GetBoundaryDensity( Flag ) * std::pow( GetBoundaryFlowRate( Flag ) / GetBoundaryArea( Flag ), 2 );
}

Real
MS_Model_1D::GetBoundaryStress( const BCFlag& Flag, const stressTypes& StressType ) const
{
    switch ( StressType )
    {
    case StaticPressure:
    {
        return GetBoundaryPressure( Flag );
    }

    case TotalPressure:
    {
        return GetBoundaryPressure( Flag ) + GetBoundaryDynamicPressure( Flag ) * ( ( GetBoundaryFlowRate( Flag ) > 0.0 ) ? 1 : -1 );
    }

    default:
    {
        std::cout << "ERROR: Invalid stress type [" << Enum2String( StressType, MS_stressesMap ) << "]" << std::endl;

        return 0.0;
    }
    }
}

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE

Real
MS_Model_1D::GetBoundaryDeltaFlowRate( const BCFlag& Flag, bool& SolveLinearSystem )
{
    OneD_BCSide bcSide = FlagConverter( Flag );

    SolveLinearModel( SolveLinearSystem );

    Real Q      = M_Solver->BoundaryValue( *M_Solution, OneD_Q, bcSide );
    Real Qdelta = M_Solver->BoundaryValue( *M_LinearSolution, OneD_Q, bcSide );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::GetBoundaryDeltaFlowRate( Flag, SolveLinearSystem ) \n";
    Debug( 8130 ) << "Q:          " << Q << "\n";
    Debug( 8130 ) << "Qdelta:     " << Qdelta << "\n";
#endif

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA

    if ( M_BCDeltaType == OneD_A )
    {
        // dQ/dP
        return ( (Qdelta - Q) / M_BCDelta ) * M_Physics->dAdP( M_Solver->BoundaryValue( *M_Solution, OneD_P, bcSide ), 0 );
    }
    else
    {
        // dQ/dQ
        return (Qdelta - Q) / M_BCDelta;
    }

#else

    return (Qdelta - Q) / M_BCDelta;

#endif

}

Real
MS_Model_1D::GetBoundaryDeltaPressure( const BCFlag& Flag, bool& SolveLinearSystem )
{
    OneD_BCSide bcSide = FlagConverter( Flag );

    SolveLinearModel( SolveLinearSystem );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA

    Real A      = M_Solver->BoundaryValue( *M_Solution, OneD_A, bcSide );
    Real Adelta = M_Solver->BoundaryValue( *M_LinearSolution, OneD_A, bcSide );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::GetBoundaryDeltaPressure( Flag, SolveLinearSystem ) \n";
    Debug( 8130 ) << "A:          " << A <<  "\n";
    Debug( 8130 ) << "Adelta:     " << Adelta <<  "\n";
#endif

    if ( M_BCDeltaType == OneD_A )
    {
        // dP/dP
        return ( (Adelta - A) / M_BCDelta ) * M_Physics->dPdA( M_Solver->BoundaryValue( *M_Solution, OneD_A, bcSide ), 0 )
               * M_Physics->dAdP( M_Solver->BoundaryValue( *M_Solution, OneD_P, bcSide ), 0 );
    }
    else
    {
        // dP/dQ
        return ( (Adelta - A) / M_BCDelta ) * M_Physics->dPdA( M_Solver->BoundaryValue( *M_Solution, OneD_A, bcSide ), 0 );
    }

#else

    Real P      = M_Solver->BoundaryValue( *M_Solution, OneD_P, bcSide );
    Real Pdelta = M_Solver->BoundaryValue( *M_LinearSolution, OneD_P, bcSide );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::GetBoundaryDeltaPressure( Flag, SolveLinearSystem ) \n";
    Debug( 8130 ) << "P:          " << P <<  "\n";
    Debug( 8130 ) << "Pdelta:     " << Pdelta <<  "\n";
#endif

    return (Pdelta - P) / M_BCDelta;

#endif

}

#else

Real
MS_Model_1D::GetBoundaryDeltaFlowRate( const BCFlag& Flag, bool& /*SolveLinearSystem*/ )
{
    return TangentProblem( FlagConverter( Flag ), OneD_Q );
}

Real
MS_Model_1D::GetBoundaryDeltaPressure( const BCFlag& Flag, bool& /*SolveLinearSystem*/ )
{
    return TangentProblem( FlagConverter( Flag ), OneD_P );
}

#endif

Real
MS_Model_1D::GetBoundaryDeltaDynamicPressure( const BCFlag& Flag, bool& SolveLinearSystem )
{
    // Untested
    return GetBoundaryDensity( Flag ) * GetBoundaryDeltaFlowRate( Flag, SolveLinearSystem ) * GetBoundaryFlowRate( Flag ) / ( GetBoundaryArea( Flag ) * GetBoundaryArea( Flag ) );
}

Real
MS_Model_1D::GetBoundaryDeltaStress( const BCFlag& Flag, bool& SolveLinearSystem, const stressTypes& StressType )
{
    switch ( StressType )
    {
    case StaticPressure:
    {
        return GetBoundaryDeltaPressure( Flag, SolveLinearSystem );
    }

    case TotalPressure:
    {
        return GetBoundaryDeltaPressure( Flag, SolveLinearSystem ) + GetBoundaryDeltaDynamicPressure( Flag, SolveLinearSystem ); //Verify the sign of DynamicPressure contribute!
    }

    default:
    {
        std::cout << "ERROR: Invalid stress type [" << Enum2String( StressType, MS_stressesMap ) << "]" << std::endl;

        return 0.0;
    }
    }
}

// ===================================================
// Get Methods
// ===================================================
MS_Model_1D::BC_Type&
MS_Model_1D::GetBC() const
{
    return *(M_BC->GetHandler());
}


MS_Model_1D::Data_Type&
MS_Model_1D::GetData() const
{
    return *M_Data;
}

MS_Model_1D::Physics_PtrType
MS_Model_1D::GetPhysics() const
{
    return M_Physics;
}

MS_Model_1D::Flux_PtrType
MS_Model_1D::GetFlux() const
{
    return M_Flux;
}

MS_Model_1D::Source_PtrType
MS_Model_1D::GetSource() const
{
    return M_Source;
}

MS_Model_1D::FESpace_PtrType
MS_Model_1D::GetFESpace() const
{
    return M_FESpace;
}

MS_Model_1D::Solver_PtrType
MS_Model_1D::GetSolver() const
{
    return M_Solver;
}

const MS_Model_1D::Solution_PtrType&
MS_Model_1D::GetSolution() const
{
    return M_Solution;
}

const MS_Model_1D::Vector_PtrType&
MS_Model_1D::GetSolution( const std::string& quantity) const
{
    return (*M_Solution)[quantity];
}

// ===================================================
// Private Methods
// ===================================================
void
MS_Model_1D::SetupGlobalData( const std::string& FileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupGlobalData( FileName ) \n";
#endif

    GetPot DataFile( FileName );

    //Global data time
    M_Data->setDataTime( M_globalData->GetDataTime() );

    //Global physical quantities
    if ( !DataFile.checkVariable( "1D_Model/PhysicalParameters/density" ) )
        M_Data->setDensity( M_globalData->GetFluidDensity() );
    if ( !DataFile.checkVariable( "1D_Model/PhysicalParameters/viscosity" ) )
        M_Data->setViscosity( M_globalData->GetFluidViscosity() );
    if ( !DataFile.checkVariable( "1D_Model/PhysicalParameters/externalPressure" ) )
        M_Data->setExternalPressure( M_globalData->GetFluidReferencePressure() );
    if ( !DataFile.checkVariable( "1D_Model/PhysicalParameters/densityWall" ) )
        M_Data->setDensityWall( M_globalData->GetStructureDensity() );
    if ( !DataFile.checkVariable( "1D_Model/PhysicalParameters/poisson" ) )
        M_Data->setPoisson( M_globalData->GetStructurePoissonCoefficient() );
    if ( !DataFile.checkVariable( "1D_Model/PhysicalParameters/young" ) )
        M_Data->setYoung( M_globalData->GetStructureYoungModulus() );

    //After changing some parameters we need to update the coefficients
    M_Data->UpdateCoefficients();
}

void
MS_Model_1D::SetupFESpace()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupFEspace() \n";
#endif

    //Transform mesh
    boost::array< Real, NDIM > NullTransformation;
    NullTransformation[0] = 0.;
    NullTransformation[1] = 0.;
    NullTransformation[2] = 0.;

    //The real mesh can be only scaled due to OneDimensionalModel_Solver conventions
    M_Data->mesh()->transformMesh( M_geometryScale, NullTransformation, NullTransformation ); // Scale the x dimension

    for ( UInt i(0); i < M_Data->NumberOfNodes() ; ++i )
        M_Data->setArea0( M_Data->Area0( i ) * M_geometryScale[1] * M_geometryScale[2], i );  // Scale the area (y-z dimensions)

    //After changing some parameters we need to update the coefficients
    M_Data->UpdateCoefficients();

#ifdef HAVE_HDF5
    //The mesh for the post-processing can be rotated
    M_ExporterMesh->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );
#endif

    //Setup FESpace
    const RefFE*    refFE = &feSegP1;
    const QuadRule* qR    = &quadRuleSeg3pt;
    const QuadRule* bdQr  = &quadRuleSeg1pt;

//    const RefFE*    refFE = &feSegP2;
//    const QuadRule* qR    = &quadRuleSeg3pt;
//    const QuadRule* bdQr  = &quadRuleSeg1pt;

    M_FESpace.reset( new FESpace_Type( M_Data->mesh(), *refFE, *qR, *bdQr, 1, M_comm ) );
    M_Solver->setFESpace( M_FESpace );
}

void
MS_Model_1D::InitializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::InitializeSolution() \n";
#endif

    if ( MS_ProblemStep > 0 )
    {
        M_Importer->setMeshProcId( M_ExporterMesh, M_comm->MyPID() );

        //M_Exporter->addVariable( ExporterData::Scalar, "Solid Area",      (*M_Solution)["A"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        M_Importer->addVariable( ExporterData::Scalar, "Area ratio",      (*M_Solution)["A/A0-1"], static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        M_Importer->addVariable( ExporterData::Scalar, "Fluid Flow Rate", (*M_Solution)["Q"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        //M_Importer->addVariable( ExporterData::Scalar, "W1",              (*M_Solution)["W1"],   static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        //M_Importer->addVariable( ExporterData::Scalar, "W2",              (*M_Solution)["W2"],   static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
        M_Importer->addVariable( ExporterData::Scalar, "Fluid Pressure",  (*M_Solution)["P"],    static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );

        // Import
        M_Exporter->setStartIndex( M_Importer->importFromTime( M_Data->dataTime()->getInitialTime() ) + 1 );

        // Compute A from AreaRatio
        M_Solver->computeArea( *M_Solution );

        // Compute W1 and W2 from A and Q
        M_Solver->computeW1W2( *M_Solution );
    }
    else
        M_Solver->initialize( *M_Solution );
}

void
MS_Model_1D::UpdateSolution( const Solution_Type& solution1, Solution_Type& solution2 )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::UpdateSolution( solution1, solution2 ) \n";
#endif

    // Here we make a true copy (not a copy of the shared_ptr, but a copy of its content)
    for ( Solution_ConstIterator i = solution1.begin() ; i != solution1.end() ; ++i )
        *solution2[i->first] = *i->second;
}

void
MS_Model_1D::Solve( BC_Type& bc, Solution_Type& solution, const std::string& solverType )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::Solve() \n";
#endif

    // Re-initialize solution
    UpdateSolution( *M_Solution_tn, solution );

    // Subiterate to respect CFL
    UInt SubiterationNumber(1);
    Real timeStep = M_Data->dataTime()->getTimeStep();

    Real CFL = M_Solver->ComputeCFL( solution, M_Data->dataTime()->getTimeStep() );
    if ( CFL > M_Data->CFLmax() )
    {
        SubiterationNumber = std::ceil( CFL / M_Data->CFLmax() );
        timeStep /= SubiterationNumber;
    }

    if ( M_displayer->isLeader() )
        std::cout << solverType << "  CFL                                      " << CFL*timeStep/M_Data->dataTime()->getTimeStep() << std::endl;

    for ( UInt i(1) ; i <= SubiterationNumber ; ++i )
    {
        if ( M_displayer->isLeader() )
        {
            std::cout << solverType << "  Subiteration                             " << i << "/" << SubiterationNumber << std::endl;
            std::cout << solverType << "  Time                                     " <<  M_Data->dataTime()->getPreviousTime() + i*timeStep << std::endl;
        }
        //bc.UpdateOperatorVariables();

        M_Solver->updateRHS( solution, timeStep );
        M_Solver->iterate( bc, solution, M_Data->dataTime()->getPreviousTime() + i*timeStep, timeStep );
    }
}

OneD_BCSide
MS_Model_1D::FlagConverter( const BCFlag& flag ) const
{
    return (flag == 0) ? OneD_left : OneD_right;
}

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE

void
MS_Model_1D::CreateLinearBC()
{
    // Allocating the correct space
    M_BCPreviousTimeSteps.reserve( std::max( M_couplings[0]->GetTimeInterpolationOrder(), M_couplings[1]->GetTimeInterpolationOrder() ) );

    // Create bcSide map
    std::map< OneD_BCSide, std::map< OneD_BC, Real > > bcSideMap;
    M_BCPreviousTimeSteps.push_back( bcSideMap );

    // Create bcType map
    std::map< OneD_BC, Real > bcTypeMap;
    M_BCPreviousTimeSteps[0][OneD_left]  = bcTypeMap;
    M_BCPreviousTimeSteps[0][OneD_right] = bcTypeMap;
}

void
MS_Model_1D::UpdateLinearBC( const Solution_Type& solution )
{
    M_BCPreviousTimeSteps[0][OneD_left][OneD_A]  = M_Solver->BoundaryValue( solution, OneD_A, OneD_left );
    M_BCPreviousTimeSteps[0][OneD_left][OneD_P]  = M_Solver->BoundaryValue( solution, OneD_P, OneD_left );
    M_BCPreviousTimeSteps[0][OneD_left][OneD_Q]  = M_Solver->BoundaryValue( solution, OneD_Q, OneD_left );
    M_BCPreviousTimeSteps[0][OneD_right][OneD_A] = M_Solver->BoundaryValue( solution, OneD_A, OneD_right );
    M_BCPreviousTimeSteps[0][OneD_right][OneD_P] = M_Solver->BoundaryValue( solution, OneD_P, OneD_right );
    M_BCPreviousTimeSteps[0][OneD_right][OneD_Q] = M_Solver->BoundaryValue( solution, OneD_Q, OneD_right );
}

void
MS_Model_1D::ImposePerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::Perturbation() \n";
#endif

    for ( MS_CouplingsVector_ConstIterator i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->IsPerturbed() )
        {
            // Find the side to perturb and apply the perturbation
            M_BCDeltaSide = FlagConverter( ( *i )->GetFlag( ( *i )->GetModelLocalID( M_ID ) ) );
            M_LinearBC->BC( M_BCDeltaSide )->setBCFunction( OneD_first, M_BCBaseDelta );

            // Compute the range
            M_BCDeltaType = M_LinearBC->BC( M_BCDeltaSide )->type( OneD_first );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA
            // We replace pressure BC with area BC for the perturbed problem
            if ( M_BCDeltaType == OneD_P )
            {
                M_LinearBC->BC( M_BCDeltaSide )->setType( OneD_first, OneD_A );
                M_BCDeltaType = OneD_A;
            }
#endif

            //M_BCDelta = ( *i )->GetResidual()[( *i )->GetPerturbedCoupling()] * 10000;
            //M_BCDelta = ( M_BCDelta[1] - M_BCDelta[0] ) / 100;

            //if ( std::abs( M_BCDelta ) < 1e-6 || std::abs( M_BCDelta ) > 1e6 )
            switch ( M_BCDeltaType )
            {
            case OneD_A:

                M_BCDelta = M_Data->JacobianPerturbationArea();
                break;

            case OneD_Q:

                M_BCDelta = M_Data->JacobianPerturbationFlowRate();

                break;

            case OneD_P:

                M_BCDelta = 5; //M_Data->JacobianPerturbationPressure();

                break;

            default:
                std::cout << "Warning: bcType \"" << M_BCDeltaType << "\"not available!" << std::endl;
            }

            break;
        }

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "BCDelta:    " << M_BCDelta << "\n";
#endif

}

void
MS_Model_1D::ResetPerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::ResetPerturbation() \n";
#endif

    M_LinearBC->BC( M_BCDeltaSide )->setBCFunction( OneD_first, M_BC->GetHandler()->BC( M_BCDeltaSide )->BCFunction( OneD_first ) );

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE_AREA
    // Restoring the original BC
    if ( M_BCDeltaType == OneD_A )
        M_LinearBC->BC( M_BCDeltaSide )->setType( OneD_first, OneD_P );
#endif

}

Real
MS_Model_1D::BCFunctionDelta( const Real& t )
{
    // Lagrange interpolation
    Real bcValue(0);
    Real base(1);

    // Time container for interpolation
    std::vector< Real > timeContainer( M_BCPreviousTimeSteps.size(), 0 );
    for ( UInt i(0) ; i < M_BCPreviousTimeSteps.size() ; ++i )
        timeContainer[i] = M_globalData->GetDataTime()->getTime() - i * M_globalData->GetDataTime()->getTimeStep();

    for ( UInt i(0) ; i < M_BCPreviousTimeSteps.size() ; ++i )
    {
        base = 1;
        for ( UInt j(0) ; j < M_BCPreviousTimeSteps.size() ; ++j )
            if ( j != i )
                base *= (t - timeContainer[j]) / (timeContainer[i] - timeContainer[j]);

        if ( i == 0 )
            bcValue += ( M_BCDelta + M_BCPreviousTimeSteps[i][M_BCDeltaSide][M_BCDeltaType] ) * base;
        else
            bcValue += M_BCPreviousTimeSteps[i][M_BCDeltaSide][M_BCDeltaType] * base;
    }

    return bcValue;
}

#else

Real
MS_Model_1D::TangentProblem( const OneD_BCSide& bcOutputSide, const OneD_BC& bcOutputType )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::TangentProblem( bcOutputSide, bcOutputType ) \n";
#endif

    Real JacobianCoefficient(0);

    for ( MS_CouplingsVector_ConstIterator i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->IsPerturbed() )
        {
            // Find the perturbed side
            OneD_BCSide bcSide = FlagConverter( ( *i )->GetFlag( ( *i )->GetModelLocalID( M_ID ) ) );

            // Perturbation has no effect on the other sides (which also means that dQ/dQ and dP/dP are always zero)
            if ( bcSide != bcOutputSide )
                break;

            // Compute the eigenvectors
            Container2D_Type eigenvalues, leftEigenvector1, leftEigenvector2;
            M_Solver->BoundaryEigenValuesEigenVectors( bcSide, *M_Solution_tn, eigenvalues, leftEigenvector1, leftEigenvector2 );

            switch ( bcSide )
            {
            case OneD_left:
                switch ( bcOutputType )
                {
                case OneD_Q: // dQ_L/dP_L
                    JacobianCoefficient = leftEigenvector2[0] / leftEigenvector2[1]
                                          * M_Physics->dAdP( M_Solver->BoundaryValue( *M_Solution, OneD_P, OneD_left ), 0 );
                    break;
                case OneD_P: // dP_L/dQ_L
                    JacobianCoefficient = leftEigenvector2[1] / leftEigenvector2[0]
                                          * M_Physics->dPdA( M_Solver->BoundaryValue( *M_Solution, OneD_A, OneD_left ), 0 );
                    break;
                default:
                    std::cout << "Warning: bcType \"" << bcOutputType << "\"not available!" << std::endl;
                }
                break;
            case OneD_right:
                switch ( bcOutputType )
                {
                case OneD_Q: // dQ_R/dP_R
                    JacobianCoefficient = -leftEigenvector1[0] / leftEigenvector1[1]
                                          * M_Physics->dAdP( M_Solver->BoundaryValue( *M_Solution, OneD_P, OneD_right ), M_Data->NumberOfElements() );
                    break;
                case OneD_P: // dP_R/dQ_R
                    JacobianCoefficient = -leftEigenvector1[1] / leftEigenvector1[0]
                                          * M_Physics->dPdA( M_Solver->BoundaryValue( *M_Solution, OneD_A, OneD_right ), M_Data->NumberOfElements() );
                    break;
                default:
                    std::cout << "Warning: bcType \"" << bcOutputType << "\"not available!" << std::endl;
                }
                break;
            default:
                std::cout << "Warning: bcSide \"" << bcSide << "\" not available!" << std::endl;
            }

            break;
        }

    return JacobianCoefficient;
}

#endif

} // Namespace LifeV
