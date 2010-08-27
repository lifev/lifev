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
 */

#include "MS_Model_1D.hpp"

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Model_1D::MS_Model_1D() :
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
    M_LinearBC                     (),
    M_Solution_tn                  ( new Solution_Type() ),
    M_Solution                     ( new Solution_Type() ),
    M_LinearSolution               ( new Solution_Type() ),
    M_FESpace                      (),
    M_LinearSolver                 (),
    M_BCBaseDelta                  (),
    M_BCDelta                      ()
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
    M_ExporterMesh                 ( OneDimensionalModel.M_ExporterMesh ),
#endif
    M_Data                         ( OneDimensionalModel.M_Data ),
    M_Physics                      ( OneDimensionalModel.M_Physics ),
    M_Flux                         ( OneDimensionalModel.M_Flux ),
    M_Source                       ( OneDimensionalModel.M_Source ),
    M_Solver                       ( OneDimensionalModel.M_Solver ),
    M_BC                           ( OneDimensionalModel.M_BC ),
    M_LinearBC                     ( OneDimensionalModel.M_LinearBC ),
    M_Solution_tn                  ( OneDimensionalModel.M_Solution_tn ),
    M_Solution                     ( OneDimensionalModel.M_Solution ),
    M_LinearSolution               ( OneDimensionalModel.M_LinearSolution ),
    M_FESpace                      ( OneDimensionalModel.M_FESpace ),
    M_LinearSolver                 ( OneDimensionalModel.M_LinearSolver ),
    M_BCBaseDelta                  ( OneDimensionalModel.M_BCBaseDelta ),
    M_BCDelta                      ( OneDimensionalModel.M_BCDelta )
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
        M_ExporterMesh                 = OneDimensionalModel.M_ExporterMesh;
#endif
        M_Data                         = OneDimensionalModel.M_Data;
        M_Physics                      = OneDimensionalModel.M_Physics;
        M_Flux                         = OneDimensionalModel.M_Flux;
        M_Source                       = OneDimensionalModel.M_Source;
        M_Solver                       = OneDimensionalModel.M_Solver;
        M_BC                           = OneDimensionalModel.M_BC;
        M_LinearBC                     = OneDimensionalModel.M_LinearBC;
        M_Solution_tn                  = OneDimensionalModel.M_Solution_tn;
        M_Solution                     = OneDimensionalModel.M_Solution;
        M_LinearSolution               = OneDimensionalModel.M_LinearSolution;
        M_FESpace                      = OneDimensionalModel.M_FESpace;
        M_LinearSolver                 = OneDimensionalModel.M_LinearSolver;
        M_BCBaseDelta                  = OneDimensionalModel.M_BCBaseDelta;
        M_BCDelta                      = OneDimensionalModel.M_BCDelta;
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
    M_BC->SetOperator( M_Solver );
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
    M_BC->GetHandler()->setDefaultBC( M_Flux, M_Source );
    M_BC->SetSolution( M_Solution );

#ifdef PERTURBATION
    SetupLinearModel();
#endif

#ifdef HAVE_HDF5
    //Post-processing
    M_Exporter->setMeshProcId( M_ExporterMesh, M_comm->MyPID() );

    M_Exporter->addVariable( ExporterData::Scalar, "Solid Area",      (*M_Solution)["A"],  static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_Exporter->addVariable( ExporterData::Scalar, "Fluid Flow Rate", (*M_Solution)["Q"],  static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_Exporter->addVariable( ExporterData::Scalar, "W1",              (*M_Solution)["W1"], static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_Exporter->addVariable( ExporterData::Scalar, "W2",              (*M_Solution)["W2"], static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
    M_Exporter->addVariable( ExporterData::Scalar, "Fluid Pressure",  (*M_Solution)["P"],  static_cast <UInt> ( 0 ), M_FESpace->dof().numTotalDof() );
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

    M_Solver->initialize( *M_Solution );
    UpdateSolution( *M_Solution, *M_Solution_tn );
}

void
MS_Model_1D::UpdateSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::UpdateSystem() \n";
#endif

    // Update previous solution
    UpdateSolution( *M_Solution, *M_Solution_tn );
}

void
MS_Model_1D::SolveSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SolveSystem() \n";
#endif

    Solve( *M_BC->GetHandler(), *M_Solution );
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

    //Matlab post-processing
    M_Solver->postProcess( *M_Solution, M_Data->dataTime()->getTime() );
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
void
MS_Model_1D::SetupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::SetupLinearModel( ) \n";
#endif

    // Solution for the linear problem (this does not change anything in the solver)
    M_Solver->setupSolution( *M_LinearSolution );

    // Define BCFunctions for tangent problem
    M_BCBaseDelta.setFunction( boost::bind( &MS_Model_1D::BCFunctionDelta, this, _1 ) );

    // The linear BCHandler is a copy of the original BCHandler with the LinearSolution instead of the true solution
    M_LinearBC.reset( new BC_Type( *M_BC->GetHandler() ) );
    M_LinearBC->setSolution( M_LinearSolution );
}

void
MS_Model_1D::UpdateLinearModel()
{
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
            return -GetBoundaryPressure( Flag );
        }

        case TotalPressure:
        {
            return -GetBoundaryPressure( Flag ) + GetBoundaryDynamicPressure( Flag ) * ( ( GetBoundaryFlowRate( Flag ) > 0.0 ) ? 1 : -1 );
        }

        default:
        {
            std::cout << "ERROR: Invalid stress type [" << Enum2String( StressType, MS_stressesMap ) << "]" << std::endl;

            return 0.0;
        }
    }
}

Real
MS_Model_1D::GetBoundaryDeltaArea( const BCFlag& Flag, bool& SolveLinearSystem )
{

#ifdef PERTURBATION
    Real A = M_Solver->BoundaryValue( *M_Solution, OneD_A, FlagConverter( Flag ) );

    SolveLinearModel( SolveLinearSystem );

    Real Adelta = M_Solver->BoundaryValue( *M_LinearSolution, OneD_A, FlagConverter( Flag ) );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::GetBoundaryDeltaArea( Flag, SolveLinearSystem ) \n";
    Debug( 8130 ) << "A:          " << A <<  "\n";
    Debug( 8130 ) << "Adelta:     " << Adelta <<  "\n";
#endif

    return (Adelta - A) / M_BCDelta[2];
#else
    return TangentProblem( FlagConverter( Flag ), OneD_A );
#endif

}

Real
MS_Model_1D::GetBoundaryDeltaFlowRate( const BCFlag& Flag, bool& SolveLinearSystem )
{

#ifdef PERTURBATION
    Real Q = M_Solver->BoundaryValue( *M_Solution, OneD_Q, FlagConverter( Flag ) );

    SolveLinearModel( SolveLinearSystem );

    Real Qdelta = M_Solver->BoundaryValue( *M_LinearSolution, OneD_Q, FlagConverter( Flag ) );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::GetBoundaryDeltaFlowRate( Flag, SolveLinearSystem ) \n";
    Debug( 8130 ) << "Q:          " << Q <<  "\n";
    Debug( 8130 ) << "Qdelta:     " << Qdelta <<  "\n";
#endif

    return (Qdelta - Q) / M_BCDelta[2];
#else
    return TangentProblem( FlagConverter( Flag ), OneD_Q );
#endif

}

Real
MS_Model_1D::GetBoundaryDeltaPressure( const BCFlag& Flag, bool& SolveLinearSystem )
{

#ifdef PERTURBATION
    Real P = M_Solver->BoundaryValue( *M_Solution, OneD_P, FlagConverter( Flag ) );

    SolveLinearModel( SolveLinearSystem );

    Real Pdelta = M_Solver->BoundaryValue( *M_LinearSolution, OneD_P, FlagConverter( Flag ) );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::GetBoundaryDeltaPressure( Flag, SolveLinearSystem ) \n";
    Debug( 8130 ) << "P:          " << P <<  "\n";
    Debug( 8130 ) << "Pdelta:     " << Pdelta <<  "\n";
#endif

    return (Pdelta - P) / M_BCDelta[2];
#else
    return TangentProblem( FlagConverter( Flag ), OneD_P );
#endif

}

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
            return -GetBoundaryDeltaPressure( Flag, SolveLinearSystem );
        }

        case TotalPressure:
        {
            return -GetBoundaryDeltaPressure( Flag, SolveLinearSystem ) + GetBoundaryDeltaDynamicPressure( Flag, SolveLinearSystem ); //Verify the sign of DynamicPressure contribute!
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
    NullTransformation[0] = 0.; NullTransformation[1] = 0.; NullTransformation[2] = 0.;

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

    M_FESpace.reset( new FESpace_Type( M_Data->mesh(), *refFE, *qR, *bdQr, 1, M_comm ) );
    M_Solver->setFESpace( M_FESpace );
}

void
MS_Model_1D::UpdateSolution( const Solution_Type& solution1, Solution_Type& solution2 )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::UpdateSolution( solution1, solution2 ) \n";
#endif

    // Here we make a true copy (not a copy of the shared_ptr, but a copy of its content)
    for ( Solution_ConstIterator i = solution1.begin() ; i != solution1.end() ; ++i )
        *solution2[i->first]= *i->second;
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

    for ( UInt i(0) ; i < SubiterationNumber ; ++i )
    {
        if ( M_displayer->isLeader() )
            std::cout << solverType << "  Subiteration                             " << (i+1) << "/" << SubiterationNumber << std::endl;

        //bc.UpdateOperatorVariables();

        M_Solver->updateRHS( solution, timeStep );
        M_Solver->iterate( bc, solution, M_Data->dataTime()->getTime() + i*timeStep, timeStep );
    }
}

OneD_BCSide
MS_Model_1D::FlagConverter( const BCFlag& flag ) const
{
    return (flag == 0) ? OneD_left : OneD_right;
}

void
MS_Model_1D::ImposePerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::ImposePerturbation() \n";
#endif

    for ( MS_CouplingsVector_ConstIterator i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->IsPerturbed() )
        {
            // Find the side to perturb and apply the perturbation
            OneD_BCSide bcFlag = FlagConverter( ( *i )->GetFlag( ( *i )->GetModelLocalID( M_ID ) ) );
            M_LinearBC->BC( bcFlag )->setBCFunction( OneD_first, M_BCBaseDelta );

            // Compute the range
            OneD_BC bcType = M_LinearBC->BC( bcFlag )->type( OneD_first );
            M_BCDelta[0] = M_Solver->BoundaryValue( *M_Solution_tn, bcType, bcFlag );
            M_BCDelta[1] = M_Solver->BoundaryValue( *M_Solution,    bcType, bcFlag );

            // Compute the perturbation (work in progress)
            M_BCDelta[2] = ( M_BCDelta[1] - M_BCDelta[0] ) * 100;
            if ( std::abs( M_BCDelta[2] ) < 1e-6)
                M_BCDelta[2] = ( *i )->GetResidual()[( *i )->GetPerturbedCoupling()];
            if ( std::abs( M_BCDelta[2] ) < 1e-6)
                M_BCDelta[2] = 1;

            std::cout << "Way1: " << ( M_BCDelta[1] - M_BCDelta[0] ) * 100 << std::endl;
            std::cout << "Way2: " << ( *i )->GetResidual()[( *i )->GetPerturbedCoupling()] << std::endl;
            std::cout << "Used: " << M_BCDelta[2] << std::endl;

            break;
        }

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "BCDelta[0]: " << M_BCDelta[0] << "\n";
    Debug( 8130 ) << "BCDelta[1]: " << M_BCDelta[1] << "\n";
    Debug( 8130 ) << "BCDelta[2]: " << M_BCDelta[2] << "\n";
#endif

}

void
MS_Model_1D::ResetPerturbation()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8130 ) << "MS_Model_1D::ResetPerturbation() \n";
#endif

    for ( MS_CouplingsVector_ConstIterator i = M_couplings.begin(); i < M_couplings.end(); ++i )
        if ( ( *i )->IsPerturbed() )
        {
            OneD_BCSide bcFlag = FlagConverter( ( *i )->GetFlag( ( *i )->GetModelLocalID( M_ID ) ) );
            M_LinearBC->BC( bcFlag )->setBCFunction( OneD_first, M_BC->GetHandler()->BC( bcFlag )->BCFunction( OneD_first ) );

            break;
        }
}

Real
MS_Model_1D::BCFunctionDelta( const Real& t )
{
    return M_BCDelta[0] + ( M_BCDelta[1] + M_BCDelta[2] - M_BCDelta[0] )
                        / ( M_Data->dataTime()->getTime() - M_Data->dataTime()->getPreviousTime() )
                        * ( t - M_Data->dataTime()->getPreviousTime() );
}

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

            // Perturbation has no effect on the other site.
            if ( bcSide != bcOutputSide )
                break;

            // Compute the eigenvectors
            Container2D_Type eigenvalues, leftEigenvector1, leftEigenvector2;
            M_Solver->BoundaryEigenValuesEigenVectors( bcSide, *M_Solution_tn, eigenvalues, leftEigenvector1, leftEigenvector2 );

            switch( bcSide )
            {
                case OneD_left:
                    switch( bcOutputType )
                    {
                        case OneD_Q:
                            JacobianCoefficient = -leftEigenvector2[0] / leftEigenvector2[1]
                                                * M_Physics->dAdP( M_Solver->BoundaryValue( *M_Solution_tn, OneD_P, bcSide ), 0 );
                            break;
                        case OneD_A:
                            JacobianCoefficient = -leftEigenvector2[1] / leftEigenvector2[0];
                            break;
                        case OneD_P:
                            JacobianCoefficient = -leftEigenvector2[1] / leftEigenvector2[0]
                                                * M_Physics->dPdA( M_Solver->BoundaryValue( *M_Solution_tn, OneD_A, bcSide ), 0 );
                            break;
                        default:
                            std::cout << "Warning: bcType \"" << bcOutputType << "\"not available!" << std::endl;
                    }
                    break;
                case OneD_right:
                    switch( bcOutputType )
                    {
                        case OneD_Q:
                            JacobianCoefficient = -leftEigenvector1[0] / leftEigenvector1[1]
                                                * M_Physics->dAdP( M_Solver->BoundaryValue( *M_Solution_tn, OneD_P, bcSide ), M_Data->NumberOfElements() );
                            break;
                        case OneD_A:
                            JacobianCoefficient = -leftEigenvector1[1] / leftEigenvector1[0];
                            break;
                        case OneD_P:
                            JacobianCoefficient = -leftEigenvector1[1] / leftEigenvector1[0]
                                                * M_Physics->dPdA( M_Solver->BoundaryValue( *M_Solution_tn, OneD_A, bcSide ), M_Data->NumberOfElements() );
                            break;
                        default:
                            std::cout << "Warning: bcType \"" << bcOutputType << "\"not available!" << std::endl;
                    }
                    break;
                default:
                    std::cout << "Warning: bcSide \"" << bcSide << "\" not available!" << std::endl;
            }

//#ifdef HAVE_LIFEV_DEBUG
            std::cout << "bcSide:              " << bcSide << "\n";
            std::cout << "bcOutputSide:        " << bcOutputSide << "\n";
            std::cout << "bcType:              " << bcOutputType << "\n";

            std::cout << "JacobianCoefficient: " << JacobianCoefficient << "\n";

            std::cout << "L11:                 " << leftEigenvector1[0] << "\n";
            std::cout << "L12:                 " << leftEigenvector1[1] << "\n";
            std::cout << "L21:                 " << leftEigenvector2[0] << "\n";
            std::cout << "L22:                 " << leftEigenvector2[1] << "\n";

            std::cout << "dAdP(1):             " << M_Physics->dAdP( M_Solver->BoundaryValue( *M_Solution_tn, OneD_P, bcSide ), 0 ) << "\n";
            std::cout << "dAdP(end):           " << M_Physics->dAdP( M_Solver->BoundaryValue( *M_Solution_tn, OneD_P, bcSide ), M_Data->NumberOfElements() ) << "\n";

            std::cout << "dPdA(1):             " << M_Physics->dPdA( M_Solver->BoundaryValue( *M_Solution_tn, OneD_A, bcSide ), 0 ) << "\n";
            std::cout << "dPdA(end):           " << M_Physics->dPdA( M_Solver->BoundaryValue( *M_Solution_tn, OneD_A, bcSide ), M_Data->NumberOfElements() ) << "\n";
//#endif

            break;
        }

    return JacobianCoefficient;
}

} // Namespace LifeV
