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
    M_exporter                     (),
    M_importer                     (),
    M_fileName                     (),
    M_fluid                        (),
    M_BC                           ( new BCInterface_Type() ),
    M_BDF                          (),
    M_data                         ( new Data_Type() ),
    M_dataMesh                     ( new DataMesh()),
    M_mesh                         (),
    M_map                          (),
    M_solution                     (),
    M_linearBC                     ( new BCInterface_Type() ),
    M_updateLinearModel            ( true ),
    M_uFESpace                     (),
    M_pFESpace                     (),
    M_lmDOF                        ( 0 ),
    M_alpha                        ( 0 ),
    M_beta                         (),
    M_RHS                          (),
    M_subiterationsMaximumNumber   (),
    M_tolerance                    (),
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
    M_fileName                     ( Fluid3D.M_fileName ),
    M_fluid                        ( Fluid3D.M_fluid ),
    M_BC                           ( Fluid3D.M_BC ),
    M_BDF                          ( Fluid3D.M_BDF ),
    M_data                         ( Fluid3D.M_data ),
    M_dataMesh                     ( Fluid3D.M_dataMesh ),
    M_mesh                         ( Fluid3D.M_mesh ),
    M_map                          ( Fluid3D.M_map ),
    M_solution                     ( Fluid3D.M_solution ),
    M_linearBC                     ( Fluid3D.M_linearBC ),
    M_updateLinearModel            ( Fluid3D.M_updateLinearModel ),
    M_uFESpace                     ( Fluid3D.M_uFESpace ),
    M_pFESpace                     ( Fluid3D.M_pFESpace ),
    M_lmDOF                        ( Fluid3D.M_lmDOF ),
    M_alpha                        ( Fluid3D.M_alpha ),
    M_beta                         ( Fluid3D.M_beta ),
    M_RHS                          ( Fluid3D.M_RHS ),
    M_subiterationsMaximumNumber   ( Fluid3D.M_subiterationsMaximumNumber ),
    M_tolerance                    ( Fluid3D.M_tolerance ),
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
        M_fileName                     = Fluid3D.M_fileName;
        M_fluid                        = Fluid3D.M_fluid;
        M_BC                           = Fluid3D.M_BC;
        M_BDF                          = Fluid3D.M_BDF;
        M_dataMesh                     = Fluid3D.M_dataMesh;
        M_mesh                         = Fluid3D.M_mesh;
        M_map                          = Fluid3D.M_map;
        M_solution                     = Fluid3D.M_solution;
        M_linearBC                     = Fluid3D.M_linearBC;
        M_updateLinearModel            = Fluid3D.M_updateLinearModel;
        M_uFESpace                     = Fluid3D.M_uFESpace;
        M_pFESpace                     = Fluid3D.M_pFESpace;
        M_lmDOF                        = Fluid3D.M_lmDOF;
        M_alpha                        = Fluid3D.M_alpha;
        M_beta                         = Fluid3D.M_beta;
        M_RHS                          = Fluid3D.M_RHS;
        M_subiterationsMaximumNumber   = Fluid3D.M_subiterationsMaximumNumber;
        M_tolerance                    = Fluid3D.M_tolerance;
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
    M_fileName = FileName;

    GetPot DataFile( FileName );

    //Fluid data
    M_data->setup( DataFile );
    if ( M_globalData.get() )
        SetupGlobalData( FileName );

    M_dataMesh->setup(DataFile, "fluid/space_discretization");

    // Parameters for the NS Iterations
    M_subiterationsMaximumNumber = DataFile( "fluid/miscellaneous/SubITMax", 0 );
    M_tolerance                  = DataFile( "fluid/miscellaneous/Tolerance", 1.e-6 );

    M_generalizedAitken.setDefaultOmega(     DataFile( "fluid/miscellaneous/Omega",        1.e-3 ) );
    M_generalizedAitken.setOmegaMin(         DataFile( "fluid/miscellaneous/range",        M_generalizedAitken.GetDefaultOmegaS()/1024, 0 ) );
    M_generalizedAitken.setOmegaMax(         DataFile( "fluid/miscellaneous/range",        M_generalizedAitken.GetDefaultOmegaS()*1024, 1 ) );
    M_generalizedAitken.UseDefaultOmega(     DataFile( "fluid/miscellaneous/fixedOmega",   false ) );
    M_generalizedAitken.setMinimizationType( DataFile( "fluid/miscellaneous/inverseOmega", true ) );

    //Boundary Conditions for the problem
    M_BC->SetOperator( M_fluid );    //MUST BE MOVED AFTER M_fluid.reset !!!
    M_BC->FillHandler( FileName, "fluid" );
    M_BC->UpdateOperatorVariables(); //MUST BE MOVED INSIDE THE UPDATE !!!

    //Setup linear problem
    SetupLinearData( FileName );

    //Setup Exporter & Importer
    SetupExporterImporter( FileName );
}

void
MS_Model_Fluid3D::SetupModel()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupProblem() \n";
#endif

    //Mesh
    SetupMesh();

    //FEspace
    SetupFEspace();

    //Add flow rate offset to BC
    M_lmDOF = M_BC->GetHandler()->getNumberBCWithType( Flux );
    SetupBCOffset( M_BC->GetHandler() );

    //Fluid
    M_fluid.reset( new Fluid_Type( *M_data, *M_uFESpace, *M_pFESpace, M_comm, M_lmDOF ) );
    GetPot DataFile( M_fileName );
    M_fluid->setUp( DataFile ); //Remove Preconditioner and Solver if possible!

    //Fluid MAP
    M_map.reset( new EpetraMap( M_fluid->getMap() ) );

    //BDF
    M_BDF.reset( new BDF_Type( M_data->dataTime()->getBDF_order() ) );

    //Problem coefficients
    M_beta.reset( new FluidVector_Type( M_map ) );
    M_RHS.reset ( new FluidVector_Type( M_map ) );

    //Post-processing
    M_exporter->setMeshProcId( M_mesh->mesh(), M_comm->MyPID() );

    M_solution.reset( new FluidVector_Type( *M_fluid->solution(), M_exporter->mapType() ) );
    if ( M_exporter->mapType() == Unique )
        M_solution->setCombineMode( Zero );

    M_exporter->addVariable( ExporterData::Vector, "Fluid Velocity", M_solution, static_cast<UInt> ( 0 ), M_uFESpace->dof().numTotalDof() );
    M_exporter->addVariable( ExporterData::Scalar, "Fluid Pressure", M_solution, 3 * M_uFESpace->dof().numTotalDof(),              M_pFESpace->dof().numTotalDof() );

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
    M_fluid->buildSystem();

    //Initialize BDF
    M_BDF->bdf_u().initialize_unk( *M_fluid->solution() );

    //Define problem coefficients
    if ( M_data->Stokes() )
    {
        M_alpha  = 0.0;
        *M_beta  = *M_fluid->solution(); //It is a stationary Navier-Stokes
        *M_RHS  *= 0.0;
    }
    else
    {
        M_alpha = M_BDF->bdf_u().coeff_der( 0 ) / M_data->dataTime()->getTimeStep();
        *M_beta = M_BDF->bdf_u().extrap();
        *M_RHS  = M_fluid->matrMass() * M_BDF->bdf_u().time_der( M_data->dataTime()->getTimeStep() );
    }

    //Set problem coefficients
    M_fluid->updateSystem( M_alpha, *M_beta, *M_RHS );
}

void
MS_Model_Fluid3D::UpdateSystem()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::UpdateSystem() \n";
#endif

    //Update BDF
    M_BDF->bdf_u().shift_right( *M_fluid->solution() );

    //Update problem coefficients
    M_alpha = M_BDF->bdf_u().coeff_der( 0 ) / M_data->dataTime()->getTimeStep();
    *M_beta = M_BDF->bdf_u().extrap();
    *M_RHS  = M_fluid->matrMass() * M_BDF->bdf_u().time_der( M_data->dataTime()->getTimeStep() );

    //Set problem coefficients
    M_fluid->updateSystem( M_alpha, *M_beta, *M_RHS );

    //Linear system need to be updated
    M_updateLinearModel = true;

    //Recompute preconditioner
    M_fluid->resetPrec( true );
}

void
MS_Model_Fluid3D::SolveSystem()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SolveSystem() \n";
#endif

    //Solve the problem
    M_fluid->iterate( *M_BC->GetHandler() );

    if ( M_subiterationsMaximumNumber > 0 )
    {
        Real residual = ( *M_beta - *M_fluid->solution() ).Norm2(); // Residual is computed on the whole solution vector;

        if ( M_displayer->isLeader() )
            std::cout << "  F-  Residual:                                " << residual << std::endl;

        M_generalizedAitken.restart();
        for ( UInt subIT = 1; subIT <= M_subiterationsMaximumNumber; ++subIT )
        {
            *M_beta += M_generalizedAitken.computeDeltaLambdaScalar( *M_beta, *M_beta - *M_fluid->solution() );

            //Linear model need to be updated!
            M_fluid->updateSystem( M_alpha, *M_beta, *M_RHS );
            M_updateLinearModel = true;

            //Solve system
            M_fluid->iterate( *M_BC->GetHandler() );

            // Check the new residual
            residual = ( *M_beta - *M_fluid->solution() ).Norm2(); // Residual is computed on the whole solution vector

            // Display subiteration information
            if ( M_displayer->isLeader() )
            {
                std::cout << "  F-  Sub-iteration n.:                        " << subIT << std::endl;
                std::cout << "  F-  Residual:                                " << residual << std::endl;
            }

            // Verify tolerance
            if ( residual <= M_tolerance )
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
    *M_solution = *M_fluid->solution();
    M_exporter->postProcess( M_data->dataTime()->getTime() );

#ifdef HAVE_HDF5
    if ( M_data->dataTime()->isLastTimeStep() )
        ( MS_DynamicCast< HDF5IOFile_Type >( M_exporter ) )->CloseFile();
#endif

}

void
MS_Model_Fluid3D::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "Velocity FE order   = " << M_data->uOrder() << std::endl
                  << "Pressure FE order   = " << M_data->pOrder() << std::endl << std::endl;

        std::cout << "Velocity DOF        = " << 3 * M_uFESpace->dof().numTotalDof() << std::endl
                  << "Pressure DOF        = " << M_pFESpace->dof().numTotalDof() << std::endl
                  << "lmDOF               = " << M_lmDOF << std::endl << std::endl;

        std::cout << "Fluid mesh maxH     = " << M_mesh->mesh()->maxH() << std::endl
                  << "Fluid mesh meanH    = " << M_mesh->mesh()->meanH() << std::endl << std::endl;

        std::cout << "NS SubITMax         = " << M_subiterationsMaximumNumber << std::endl
                  << "NS Tolerance        = " << M_tolerance << std::endl << std::endl << std::endl << std::endl;
    }
}


// ===================================================
// Methods
// ===================================================
void
MS_Model_Fluid3D::SetupLinearData( const std::string& FileName )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupLinearData( FileName ) \n";
#endif

    // Boundary Conditions for the linear problem
    M_linearBC->SetOperator( M_fluid );
    M_linearBC->FillHandler( FileName, "linear_fluid" );
}

void
MS_Model_Fluid3D::SetupLinearModel()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupLinearModel( ) \n";
#endif

    SetupBCOffset( M_linearBC->GetHandler() );
}

void
MS_Model_Fluid3D::UpdateLinearModel()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::UpdateLinearModel() \n";
#endif

    //Create an empty vector
    FluidVector_Type VectorZero( *M_solution ); VectorZero = 0.0;

    //UpdateLinearModel
    M_fluid->updateLinearSystem( M_fluid->matrNoBC(),
                                 M_alpha,
                                 *M_beta,
                                 *M_fluid->solution(),
                                 VectorZero,
                                 VectorZero,
                                 VectorZero,
                                 VectorZero );

    //Update Properties of BC
    M_linearBC->UpdateOperatorVariables();

    //Linear System Updated
    M_updateLinearModel = false;
}

void
MS_Model_Fluid3D::SolveLinearModel( bool& SolveLinearSystem )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SolveLinearModel() \n";
#endif

    if ( !SolveLinearSystem )
        return;

    if ( M_updateLinearModel )
        UpdateLinearModel();

    //Solve the linear problem
    M_fluid->iterateLin( *M_linearBC->GetHandler() );

    //This flag avoid recomputation of the same system
    SolveLinearSystem = false;
}

// ===================================================
// Set Methods
// ===================================================
void
MS_Model_Fluid3D::SetSolution( const boost::shared_ptr< FluidVector_Type >& Solution )
{
    M_solution = Solution;

    SetupSolution();
}

// ===================================================
// Get Methods (couplings)
// ===================================================
MS_Model_Fluid3D::BCInterface_Type&
MS_Model_Fluid3D::GetBCInterface()
{
    return *M_BC;
}

MS_Model_Fluid3D::BCInterface_Type&
MS_Model_Fluid3D::GetLinearBCInterface()
{
    return *M_linearBC;
}

Real
MS_Model_Fluid3D::GetBoundaryDensity( const BCFlag& /*Flag*/ ) const
{
    return M_data->density();
}

Real
MS_Model_Fluid3D::GetBoundaryViscosity( const BCFlag& /*Flag*/ ) const
{
    return M_data->viscosity();
}

Real
MS_Model_Fluid3D::GetBoundaryArea( const BCFlag& Flag ) const
{
    return M_fluid->area( Flag );
}

Real
MS_Model_Fluid3D::GetBoundaryFlowRate( const BCFlag& Flag ) const
{
    return M_fluid->flux( Flag );
}

Real
MS_Model_Fluid3D::GetBoundaryPressure( const BCFlag& Flag ) const
{
    return M_fluid->pressure( Flag );
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
    return M_fluid->LagrangeMultiplier(Flag, *M_BC->GetHandler() );
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

            std::cout << "ERROR: Invalid stress type [" << Enum2String( StressType, MS_stressesMap ) << "]" << std::endl;

            return 0.0;
    }
}

Real
MS_Model_Fluid3D::GetBoundaryDeltaFlux( const BCFlag& Flag, bool& SolveLinearSystem )
{
    SolveLinearModel( SolveLinearSystem );

    return M_fluid->GetLinearFlux( Flag );
}

Real
MS_Model_Fluid3D::GetBoundaryDeltaPressure( const BCFlag& Flag, bool& SolveLinearSystem )
{
    SolveLinearModel( SolveLinearSystem );

    return M_fluid->GetLinearPressure( Flag );
}

Real
MS_Model_Fluid3D::GetBoundaryDeltaDynamicPressure( const BCFlag& Flag, bool& SolveLinearSystem )
{
    SolveLinearModel( SolveLinearSystem );

    return GetBoundaryDensity( Flag ) * M_fluid->GetLinearFlux( Flag ) * GetBoundaryFlowRate( Flag ) / ( GetBoundaryArea( Flag ) * GetBoundaryArea( Flag ) );
}

Real
MS_Model_Fluid3D::GetBoundaryDeltaLagrangeMultiplier( const BCFlag& Flag, bool& SolveLinearSystem )
{
    SolveLinearModel( SolveLinearSystem );

    return M_fluid->LinearLagrangeMultiplier(Flag, *M_linearBC->GetHandler() );
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

            std::cout << "ERROR: Invalid stress type [" << Enum2String( StressType, MS_stressesMap ) << "]" << std::endl;

            return 0.0;
    }
}

// ===================================================
// Get Methods
// ===================================================
const MS_Model_Fluid3D::Data_Type&
MS_Model_Fluid3D::GetData() const
{
    return *M_data;
}

const MS_Model_Fluid3D::FluidVector_Type&
MS_Model_Fluid3D::GetSolution() const
{
    return *M_solution;
}

// ===================================================
// Private Methods
// ===================================================
void
MS_Model_Fluid3D::SetupGlobalData( const std::string& FileName )
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupGlobalData( FileName ) \n";
#endif

    GetPot DataFile( FileName );

    //Global data time
    M_data->setDataTime( M_globalData->GetDataTime() );

    //Global physical quantities
    if ( !DataFile.checkVariable( "fluid/physics/density" ) )
        M_data->density( M_globalData->GetFluidDensity() );
    if ( !DataFile.checkVariable( "fluid/physics/viscosity" ) )
        M_data->viscosity( M_globalData->GetFluidViscosity() );
}

void
MS_Model_Fluid3D::SetupExporterImporter( const std::string& FileName )
{
    GetPot DataFile( FileName );

    //Exporter
    const std::string exporterType = DataFile( "exporter/type", "ensight" );

#ifdef HAVE_HDF5
    if ( !exporterType.compare( "hdf5" ) )
        M_exporter.reset( new HDF5IOFile_Type() );
    else
#endif
        M_exporter.reset( new EnsightIOFile_Type() );

    M_exporter->setDataFromGetPot( DataFile );
    M_exporter->setPrefix( "Step_" + number2string( MS_ProblemStep ) + "_Model_" + number2string( M_ID ) );
    M_exporter->setDirectory( MS_ProblemFolder );

    //Importer
    const std::string importerType = DataFile( "importer/type", "ensight" );

#ifdef HAVE_HDF5
    if ( !importerType.compare( "hdf5" ) )
        M_importer.reset( new HDF5IOFile_Type() );
    else
#endif
        M_importer.reset( new EnsightIOFile_Type() );

    M_importer->setDataFromGetPot( DataFile );
    M_importer->setPrefix( "Step_" + number2string( MS_ProblemStep - 1 ) + "_Model_" + number2string( M_ID ) );
    M_importer->setDirectory( MS_ProblemFolder );
}

void
MS_Model_Fluid3D::SetupMesh()
{
    //Read fluid mesh from file
    boost::shared_ptr<Mesh_Type> fluidMesh(new Mesh_Type);
    readMesh( *fluidMesh, *M_dataMesh);

    //Transform mesh
    fluidMesh->transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    //Partition mesh
    M_mesh.reset( new PartitionMesh_Type( fluidMesh, M_comm ) );
}

void
MS_Model_Fluid3D::SetupFEspace()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupFEspace() \n";
#endif

    //Velocity FE Space
    const RefFE* u_refFE;
    const QuadRule* u_qR;
    const QuadRule* u_bdQr;

    if ( M_data->uOrder().compare( "P2" ) == 0 )
    {
        u_refFE = &feTetraP2;
        u_qR = &quadRuleTetra15pt; // DoE 5
        u_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else
        if ( M_data->uOrder().compare( "P1" ) == 0 )
        {
            u_refFE = &feTetraP1;
            u_qR = &quadRuleTetra4pt; // DoE 2
            u_bdQr = &quadRuleTria3pt; // DoE 2
        }
        else
            if ( M_data->uOrder().compare( "P1Bubble" ) == 0 )
            {
                u_refFE = &feTetraP1bubble;
                u_qR = &quadRuleTetra64pt; // DoE 2
                u_bdQr = &quadRuleTria3pt; // DoE 2
            }
            else
            {
                if ( M_displayer->isLeader() )
                    std::cout << M_data->uOrder() << " Velocity FE not implemented yet." << std::endl;
                exit( EXIT_FAILURE );
            }

    //Pressure FE Space
    const RefFE* p_refFE;
    const QuadRule* p_qR;
    const QuadRule* p_bdQr;

    if ( M_data->pOrder().compare( "P2" ) == 0 )
    {
        p_refFE = &feTetraP2;
        p_qR = u_qR;
        p_bdQr = &quadRuleTria3pt; // DoE 2
    }
    else
        if ( M_data->pOrder().compare( "P1" ) == 0 )
        {
            p_refFE = &feTetraP1;
            p_qR = u_qR;
            p_bdQr = &quadRuleTria3pt; // DoE 2
        }
        else
        {
            if ( M_displayer->isLeader() )
                std::cout << M_data->pOrder() << " pressure FE not implemented yet." << std::endl;
            exit( EXIT_FAILURE );
        }

    M_uFESpace.reset( new FESpace_Type( *M_mesh, *u_refFE, *u_qR, *u_bdQr, 3, M_comm ) );
    M_pFESpace.reset( new FESpace_Type( *M_mesh, *p_refFE, *p_qR, *p_bdQr, 1, M_comm ) );
}

void
MS_Model_Fluid3D::SetupDOF()
{

#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::SetupDOF \n";
#endif

    M_lmDOF = M_BC->GetHandler()->getNumberBCWithType( Flux );

    //M_uDOF = M_uFESpace->map().getMap(Unique)->NumGlobalElements();
    //M_pFESpace->dof().numTotalDof() = M_pFESpace->map().getMap(Unique)->NumGlobalElements();
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
        M_importer->setMeshProcId( M_mesh->mesh(), M_comm->MyPID() );

        M_importer->addVariable( ExporterData::Vector, "Velocity", M_solution, static_cast <UInt> ( 0 ),            M_uFESpace->dof().numTotalDof() );
        M_importer->addVariable( ExporterData::Scalar, "Pressure", M_solution, 3 * M_uFESpace->dof().numTotalDof(), M_pFESpace->dof().numTotalDof());

        // Import
        M_exporter->setStartIndex( M_importer->importFromTime( M_data->dataTime()->getInitialTime() ) + 1 );
    }
    else
        *M_solution = 0.0;

    M_fluid->initialize( *M_solution );
}

} // Namespace LifeV
