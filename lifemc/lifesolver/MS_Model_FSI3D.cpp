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
 *  @brief MultiScale Model FSI3D
 *
 *  @author Paolo Crosetto <paolo.crosetto@epfl.ch>
 *  @date 19-04-2010
 */

#include <lifemc/lifesolver/MS_Model_FSI3D.hpp>
#include <lifemc/lifesolver/Monolithic.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Model_FSI3D::MS_Model_FSI3D() :
    super                          (),
    M_solver                       (),
    M_data                         ( new data_Type() ),
    M_exporterFluid                (),
    M_exporterSolid                (),
    M_fluidVelocityPressure        (),
    M_fluidDisplacement            (),
    M_solidVelocity                (),
    M_solidDisplacement            (),
    M_fluidBC                      ( new BCInterface_Type() ),
    M_solidBC                      ( new BCInterface_Type() ),
    M_harmonicExtensionBC          ( new BCInterface_Type() ),
    M_linearizedFluidBC            (),
    M_linearizedSolidBC            ()
{

#ifdef DEBUG
    Debug( 8140 ) << "MS_Model_FSI3D::MS_Model_FSI3D() \n";
#endif

    M_type = FSI3D;
}

MS_Model_FSI3D::~MS_Model_FSI3D()
{

}

// ===================================================
// Methods
// ===================================================
void
MS_Model_FSI3D::SetupData( const std::string& fileName )
{
    super::SetupData( fileName );

    GetPot dataFile( fileName );

    // Setup data
    M_data->setup( dataFile );

    if ( M_globalData.get() )
        SetupGlobalData( fileName );

    M_solver.reset( new FSISolver_Type( M_data->method() ) );
    M_solver->setData( M_data );
    M_solver->FSIOper()->setDataFile( dataFile );

    // Mesh transformation
    M_solver->FSIOper()->fluidMesh().transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );
    M_solver->FSIOper()->solidMesh().transformMesh( M_geometryScale, M_geometryRotate, M_geometryTranslate );

    // ##################################################
    // ### THIS BLOCK SHOULD BE MOVED IN SetupModel() ###
    // ##################################################
    // Setup FEspace & DOF
    M_solver->FSIOper()->setupFEspace();
    M_solver->FSIOper()->setupDOF();
    M_solver->FSIOper()->setupFluidSolid();
    // ##################################################

    // Setup Boundary Conditions
	setupBC( fileName );

	// Setup Boundary Conditions for segregated FSI
    if( !M_data->isMonolithic() )
        setupSegregatedBC( fileName );

	// Setup exporter
    if( M_solver->isFluid() )
        SetupExporter( M_exporterFluid, dataFile );
    if( M_solver->isSolid() )
        SetupExporter( M_exporterSolid, dataFile, "_Solid" );
}

void
MS_Model_FSI3D::SetupModel()
{
    if( M_solver->isFluid() )
        SetExporterFluid( M_exporterFluid );

    if( M_solver->isSolid() )
        SetExporterSolid( M_exporterSolid );

    //to kill when everithing is merged
    //M_offset=(dynamic_cast<LifeV::Monolithic*>(M_solver->FSIOper().get()))->getOffset();
}

void
MS_Model_FSI3D::BuildSystem()
{
    M_solver->FSIOper()->setupSystem();
    M_solver->FSIOper()->buildSystem();
}

void
MS_Model_FSI3D::UpdateSystem()
{
    // Update BC
    if( M_solver->isFluid() )
    {
        M_fluidBC->UpdateOperatorVariables();
        M_harmonicExtensionBC->UpdateOperatorVariables();

        if( !M_data->isMonolithic() )
            M_linearizedFluidBC->UpdateOperatorVariables();
    }
    else
    {
        M_solidBC->UpdateOperatorVariables();

        if( !M_data->isMonolithic() )
            M_linearizedSolidBC->UpdateOperatorVariables();
    }
}

void
MS_Model_FSI3D::SolveSystem()
{
    if( M_solver->isSolid() )
    {
        M_solver->FSIOper()->getSolidDisp( *M_solidDisplacement );//    displacement(), M_offset);
        //*M_solidDisplacement *= M_solver->timeStep()*M_solver->FSIOper()->solid().rescaleFactor();
        M_solver->FSIOper()->getSolidVel( *M_solidVelocity );//    displacement(), M_offset);
    }

    if( M_solver->isFluid() )
        M_solver->FSIOper()->getFluidVelAndPres( *M_fluidVelocityPressure );

    //Solve the problem
    M_solver->iterate();

    if( M_solver->isFluid() )
        *M_fluidDisplacement = M_solver->FSIOper()->meshDisp();
    //M_time += M_solver->timeStep();
    // NOTE: the last iteration will not be post-processed in this way
}

void
MS_Model_FSI3D::SaveSolution()
{
    // Post-process must be made here. We need to add to HDF5exporter some methods to:
    // 1) save only the xmf file (done once, as it is independent from the variable)
    // 2) save a specific variable at a specific time (defined by variablename && time)
/*
    // save xmf file for the fluid at the  present time step
    if( M_solver->isFluid() )
    {
        // if ( segregated && semiimplicit false ) || fullMonolithic
        // save fluidDisplacement at this time step

        // if ( segregated && semiimplicit true)  || monolithic
        // save fluidDisplacement at the previous time step

        // save velocity and pressure always at this time step
    }

    // For the solid we can use classical approach
    if( M_solver->isSolid() )
    {
        M_solver->FSIOper()->getSolidDisp( *M_solidDisplacement );// displacement(), M_offset);
        M_solver->FSIOper()->getSolidVel( *M_solidVelocity );// displacement(), M_offset);
    }
*/

    if( M_solver->isFluid() )
        M_exporterFluid->postProcess( M_data->dataFluid()->dataTime()->getTime() );
    if( M_solver->isSolid() )
        M_exporterSolid->postProcess( M_data->dataSolid()->dataTime()->getTime() );

#ifdef HAVE_HDF5
     if ( M_data->dataFluid()->dataTime()->isLastTimeStep() )
     {
         if( M_solver->isFluid() )
             ( MS_DynamicCast< hdf5IOFile_Type >( M_exporterFluid ) )->CloseFile();
         if( M_solver->isSolid() )
             ( MS_DynamicCast< hdf5IOFile_Type >( M_exporterSolid ) )->CloseFile();
     }
#endif

}

void
MS_Model_FSI3D::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "FSI method          = " << M_data->method() << std::endl << std::endl;

        std::cout << "Velocity FE order   = " << M_solver->FSIOper()->dataFluid().uOrder() << std::endl
                  << "Pressure FE order   = " << M_solver->FSIOper()->dataFluid().pOrder() << std::endl
                  << "Structure FE order  = " << M_solver->FSIOper()->dataSolid().order() << std::endl<< std::endl;

        std::cout << "Velocity DOF        = " << 3 * M_solver->FSIOper()->uFESpace().dof().numTotalDof() << std::endl
                  << "Pressure DOF        = " << M_solver->FSIOper()->pFESpace().dof().numTotalDof() << std::endl
                  //<< "Lagrange m. DOF     = " << M_lmDOF << std::endl
                  << "Harmonic ext. DOF   = " << M_solver->FSIOper()->mmFESpace().dof().numTotalDof() << std::endl
                  << "Structure DOF       = " << M_solver->FSIOper()->dFESpace().dof().numTotalDof() << std::endl << std::endl;

        std::cout << "Fluid mesh maxH     = " << M_solver->FSIOper()->uFESpace().mesh()->maxH() << std::endl
                  << "Fluid mesh meanH    = " << M_solver->FSIOper()->uFESpace().mesh()->meanH() << std::endl
                  << "Solid mesh maxH     = " << M_solver->FSIOper()->dFESpace().mesh()->maxH() << std::endl
                  << "Solid mesh meanH    = " << M_solver->FSIOper()->dFESpace().mesh()->meanH() << std::endl << std::endl;
    }
}

// ===================================================
// Get Methods
// ===================================================
MS_Model_FSI3D::BCInterface_Type&
MS_Model_FSI3D::GetBCInterface()
{
    return *M_fluidBC;
}

MS_Model_FSI3D::BCInterface_Type&
MS_Model_FSI3D::GetLinearBCInterface()
{
    //return *M_LinearFluidBC;
}

Real
MS_Model_FSI3D::GetBoundaryDensity( const BCFlag& /*Flag*/ ) const
{
    return M_solver->FSIOper()->dataFluid().density();
}

Real
MS_Model_FSI3D::GetBoundaryViscosity( const BCFlag& /*Flag*/ ) const
{
    return M_solver->FSIOper()->dataFluid().viscosity();
}

Real
MS_Model_FSI3D::GetBoundaryArea( const BCFlag& Flag ) const
{
    return M_solver->FSIOper()->fluid().area( Flag );
}

Real
MS_Model_FSI3D::GetBoundaryFlowRate( const BCFlag& Flag ) const
{
    return M_solver->FSIOper()->fluid().flux( Flag, *M_solver->FSIOper()->un() );
}

Real
MS_Model_FSI3D::GetBoundaryPressure( const BCFlag& Flag ) const
{
    return M_solver->FSIOper()->fluid().pressure( Flag, *M_solver->FSIOper()->un() );
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
    //return M_Fluid->LagrangeMultiplier(Flag, *M_fluidBC->GetHandler(), un() ); // Overload of function to pass the vector of the solution
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
MS_Model_FSI3D::GetBoundaryDeltaFlux( const BCFlag& Flag, bool& SolveLinearSystem )
{
    //SolveLinearModel( SolveLinearSystem );

    //return M_Fluid->GetLinearFlux( Flag );
}

Real
MS_Model_FSI3D::GetBoundaryDeltaPressure( const BCFlag& Flag, bool& SolveLinearSystem )
{
    //SolveLinearModel( SolveLinearSystem );

    //return M_Fluid->GetLinearPressure( Flag );
}

Real
MS_Model_FSI3D::GetBoundaryDeltaDynamicPressure( const BCFlag& Flag, bool& SolveLinearSystem )
{
    //SolveLinearModel( SolveLinearSystem );

    //return GetBoundaryDensity( Flag ) * M_Fluid->GetLinearFlux( Flag ) * GetBoundaryFlowRate( Flag ) / ( GetBoundaryArea( Flag ) * GetBoundaryArea( Flag ) );
}

Real
MS_Model_FSI3D::GetBoundaryDeltaLagrangeMultiplier( const BCFlag& Flag, bool& SolveLinearSystem )
{
    //SolveLinearModel( SolveLinearSystem );

    //return M_Fluid->LinearLagrangeMultiplier(Flag, *M_LinearFluidBC->GetHandler() );
}

Real
MS_Model_FSI3D::GetBoundaryDeltaStress( const BCFlag& Flag, bool& SolveLinearSystem, const stressTypes& StressType )
{
    /*
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
    */
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
MS_Model_FSI3D::SetExporterFluid( const IOFile_PtrType& exporter )
{
    M_fluidVelocityPressure.reset( new vector_Type( M_solver->FSIOper()->fluid().getMap(),  M_exporterFluid->mapType() ) );
    M_fluidDisplacement.reset    ( new vector_Type( M_solver->FSIOper()->mmFESpace().map(), M_exporterFluid->mapType() ) );

    exporter->setMeshProcId( M_solver->FSIOper()->uFESpace().mesh(),
                             M_solver->FSIOper()->uFESpace().map().Comm().MyPID() );

    exporter->addVariable( ExporterData::Vector, "Fluid Velocity", M_fluidVelocityPressure, UInt(0), M_solver->FSIOper()->uFESpace().dof().numTotalDof()  );
    exporter->addVariable( ExporterData::Scalar, "Fluid Pressure", M_fluidVelocityPressure, UInt(3 * M_solver->FSIOper()->uFESpace().dof().numTotalDof()  ),
                                                                                            UInt(    M_solver->FSIOper()->pFESpace().dof().numTotalDof()  ) );
    exporter->addVariable( ExporterData::Vector, "Fluid Displacement", M_fluidDisplacement, UInt(0), M_solver->FSIOper()->mmFESpace().dof().numTotalDof() );
}

void
MS_Model_FSI3D::SetExporterSolid( const IOFile_PtrType& exporter )
{
    M_solidDisplacement.reset( new vector_Type( M_solver->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ) );
    M_solidVelocity.reset    ( new vector_Type( M_solver->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ) );

    exporter->setMeshProcId( M_solver->FSIOper()->dFESpace().mesh(),
                             M_solver->FSIOper()->dFESpace().map().Comm().MyPID() );

    exporter->addVariable( ExporterData::Vector, "Solid Velocity",     M_solidVelocity,     UInt(0), M_solver->FSIOper()->dFESpace().dof().numTotalDof() );
    exporter->addVariable( ExporterData::Vector, "Solid Displacement", M_solidDisplacement, UInt(0), M_solver->FSIOper()->dFESpace().dof().numTotalDof() );
}

void
MS_Model_FSI3D::setupBC( const std::string& fileName )
{
    // Setup BCHandler for FSI
    M_fluidBC->SetOperator( M_solver->FSIOper() );
    M_solidBC->SetOperator( M_solver->FSIOper() );
    M_harmonicExtensionBC->SetOperator( M_solver->FSIOper() );

    if ( M_solver->isFluid() )
    {
        M_fluidBC->FillHandler( fileName, "fluid" );
        M_harmonicExtensionBC->FillHandler( fileName, "mesh_motion" );
    }

    if ( M_solver->isSolid() )
        M_solidBC->FillHandler( fileName, "solid" );

    M_solver->setFluidBC( M_fluidBC->GetHandler() );
    M_solver->setSolidBC( M_solidBC->GetHandler() );
    M_solver->setHarmonicExtensionBC( M_harmonicExtensionBC->GetHandler() );

    /*update bc without bcInterface*/
    //M_solver->setFluxBC(BCh_monolithicFlux());
    //M_solver->setup(/*data_file*/);
    //M_solver->setRobinBC(BCh_monolithicRobin(*M_solver->FSIOper()));
    //M_solver->setFluidBC(BCh_monolithicFluid(*M_solver->FSIOper()));
    //M_solver->setHarmonicExtensionBC (BCh_harmonicExtension(*M_solver->FSIOper()));
    //M_solver->setSolidBC(BCh_monolithicSolid(*M_solver->FSIOper()));
}

void
MS_Model_FSI3D::setupSegregatedBC( const std::string& fileName )
{
    // Setup additional BCHandler for segregated FSI
    M_linearizedFluidBC.reset( new BCInterface_Type() );
    M_linearizedSolidBC.reset( new BCInterface_Type() );

    M_linearizedFluidBC->SetOperator( M_solver->FSIOper() );
    M_linearizedSolidBC->SetOperator( M_solver->FSIOper() );

    if ( M_solver->isFluid() )
        M_linearizedFluidBC->FillHandler( fileName, "lin_fluid" );

    if ( M_solver->isSolid() )
        M_linearizedSolidBC->FillHandler( fileName, "lin_solid" );

    M_solver->setLinFluidBC( M_linearizedFluidBC->GetHandler() );
    M_solver->setLinSolidBC( M_linearizedSolidBC->GetHandler() );
}

} // Namespace LifeV
