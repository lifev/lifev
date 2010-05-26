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
    @file
    @brief A short description of the file content

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 19 Apr 2010

    A more detailed description of the file (if necessary)
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
    M_exporterSolid                (),
    M_exporterFluid                (),
    M_solidDisp(),
    M_solidVel(),
    M_velAndPressure(),
    M_fluidDisp()
{
#ifdef DEBUG
    Debug( 8120 ) << "MS_Model_Fluid3D::MS_Model_Fluid3D() \n";
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
MS_Model_FSI3D::SetupData( const std::string& FileName )
{
    super::SetupData(FileName);
    GetPot dataFile(FileName);
    M_solver.reset(new FSISolverType(dataFile("problem/method", "monolithic")));
	M_solver->setDataFromGetPot(dataFile);

    if(M_solver->isFluid())
        SetUpExporter(M_exporterFluid, dataFile( "exporter/fluid/filename", "fluid"), dataFile );
    if(M_solver->isSolid())
        SetUpExporter(M_exporterSolid, dataFile( "exporter/solid/filename", "solid"), dataFile );
    //M_time = dataFile("fluid/time_discretization/initialtime"   ,0.);
    if(M_solver->isMonolithic())
        updateBCMonolithic(FileName);
    else
        updateBCSegregated(FileName);
}

void
MS_Model_FSI3D::SetupGlobalData( const boost::shared_ptr< MS_PhysicalData >& PhysicalData )
{

    //PAOLO: this new method is for setting global variables in place of
    //       standard/default values. It is mandatory for global dataTime.
    //       this method has to be called after SetupData( dataFile )
    //       In the future I will improve it. For now it is required.

    //Global data time
    //M_Data->setDataTime( PhysicalData->GetDataTime() );

    //Global physical quantities
    //M_Data->setDensity( PhysicalData->GetFluidDensity() );
    //M_Data->setViscosity( PhysicalData->GetFluidViscosity() );
    //M_Data->setPoisson( PhysicalData->GetStructurePoissonCoefficient() );
    //M_Data->setYoung( PhysicalData->GetStructureYoungModulus() );
}

void
MS_Model_FSI3D::SetupModel()
{
    /*initialization*/
    if(M_solver->isFluid())
        SetExporterFluid(M_exporterFluid);
    if(M_solver->isSolid())
        SetExporterSolid(M_exporterSolid);
    //to kill when everithing is merged
    //M_offset=(dynamic_cast<LifeV::Monolithic*>(M_solver->FSIOper().get()))->getOffset();
}

void
MS_Model_FSI3D::BuildSystem()
{
    /*already done in SetupModel*/
}

void
MS_Model_FSI3D::UpdateSystem()
{
    /*already done in SolveSystem*/
        if(M_solver->isFluid())
        {
            M_fluidBC->UpdateOperatorVariables();
            M_harmonicExtensionBC->UpdateOperatorVariables();
            if(!M_solver->isMonolithic())
                M_linearizedFluidBC->UpdateOperatorVariables();
        }
        else
        {
            M_solidBC->UpdateOperatorVariables();
            if(!M_solver->isMonolithic())
                M_linearizedSolidBC->UpdateOperatorVariables();
        }
}

void
MS_Model_FSI3D::SolveSystem()
{
    if(M_solver->isSolid())
    {
        M_solver->FSIOper()->getSolidDisp(*M_solidDisp);//    displacement(), M_offset);
        //*M_solidDisp *= M_solver->timeStep()*M_solver->FSIOper()->solid().rescaleFactor();
        M_solver->FSIOper()->getSolidVel(*M_solidVel);//    displacement(), M_offset);
    }
    if(M_solver->isFluid())
        M_solver->FSIOper()->getFluidVelAndPres(*M_velAndPressure);
    Real time = M_solver->FSIOper()->dataFluid().dataTime()->getTime()+M_solver->FSIOper()->dataFluid().dataTime()->getTimeStep();
        M_solver->iterate( time );
    if(M_solver->isFluid())
        *M_fluidDisp      = M_solver->FSIOper()->meshDisp();
    //M_time += M_solver->timeStep();
    // NOTE: the last iteration will not be post-processed in this way
}

void
MS_Model_FSI3D::SaveSolution()
{
    /*post-process*/
    if(M_solver->isSolid())
        M_exporterSolid->postProcess( M_solver->FSIOper()->dataFluid().dataTime()->getTime() );
    if(M_solver->isFluid())
        M_exporterFluid->postProcess( M_solver->FSIOper()->dataFluid().dataTime()->getTime() );

// #ifdef HAVE_HDF5
//     if ( M_solver->FSIOper()->dataFluid().dataTime()->isLastTimeStep() )
//     {
//         if(M_solver->isSolid())
//             M_exporterSolid->CloseFile();
//         if(M_solver->isFluid())
//             M_exporterFluid->CloseFile();
//     }
// #endif

}

void
MS_Model_FSI3D::ShowMe()
{
    /*show-me*/
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
    return *M_fluidBC;//    return *M_LinearFluidBC;
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
    return M_solver->FSIOper()->fluid().flux( Flag, *M_solver->FSIOper()->un() ); //VErify that is the flux
}

Real
MS_Model_FSI3D::GetBoundaryPressure( const BCFlag& Flag ) const
{
    return M_solver->FSIOper()->fluid().pressure( Flag, *M_solver->FSIOper()->un() ); //VErify that is the scalar value of the pressure
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

            std::cout << "ERROR: Invalid stress type [" << Enum2String( StressType, stressMap ) << "]" << std::endl;

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

            std::cout << "ERROR: Invalid stress type [" << Enum2String( StressType, stressMap ) << "]" << std::endl;

            return 0.0;
    }
    */
}

// ===================================================
// Private Methods
// ===================================================

void
MS_Model_FSI3D::SetUpExporter( filter_ptrtype& exporter, const std::string name, const GetPot& dataFile )
{
    std::string const exporterType =  dataFile( "exporter/type", "ensight");
#ifdef HAVE_HDF5
    if (exporterType.compare("hdf5") == 0)
    {
        exporter.reset( new  hdf5filter_type( dataFile, name) );
    }
    else
#endif
    {
        if (exporterType.compare("none") == 0)
        {
            exporter.reset( new NoExport<mesh_type> ( dataFile, name) );
        } else {
            exporter.reset( new  ensightfilter_type( dataFile, name) );
        }
    }
    exporter->setPrefix( "Step_" + number2string( MS_ProblemStep ) + "_Model_" + number2string( M_ID )+name );
    exporter->setDirectory( MS_ProblemFolder );

}


void
MS_Model_FSI3D::SetExporterFluid(const filter_ptrtype& exporter)
{
    M_velAndPressure.reset( new vector_type( M_solver->FSIOper()->fluid().getMap(), M_exporterFluid->mapType() ));
    M_fluidDisp.reset     ( new vector_type( M_solver->FSIOper()->mmFESpace().map(), M_exporterFluid->mapType() ));

    exporter->setMeshProcId(M_solver->FSIOper()->uFESpace().mesh(), M_solver->FSIOper()->uFESpace().map().Comm().MyPID());
    exporter->addVariable( ExporterData::Vector, "f-velocity", M_velAndPressure,
                                  UInt(0), M_solver->FSIOper()->uFESpace().dof().numTotalDof() );
    exporter->addVariable( ExporterData::Scalar, "f-pressure", M_velAndPressure,
                                  UInt(3*M_solver->FSIOper()->uFESpace().dof().numTotalDof()),
                                  UInt(M_solver->FSIOper()->pFESpace().dof().numTotalDof()) );
    exporter->addVariable( ExporterData::Vector, "f-displacement", M_fluidDisp,
                                  UInt(0), M_solver->FSIOper()->mmFESpace().dof().numTotalDof() );
}


void
MS_Model_FSI3D::SetExporterSolid(const filter_ptrtype& exporter)
{
    M_solidDisp.reset( new vector_type( M_solver->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));
    M_solidVel.reset ( new vector_type( M_solver->FSIOper()->dFESpace().map(), M_exporterSolid->mapType() ));

    exporter->setMeshProcId(M_solver->FSIOper()->dFESpace().mesh(), M_solver->FSIOper()->dFESpace().map().Comm().MyPID());

    exporter->addVariable( ExporterData::Vector, "s-displacement", M_solidDisp,
                           UInt(0), M_solver->FSIOper()->dFESpace().dof().numTotalDof() );
    exporter->addVariable( ExporterData::Vector, "s-velocity", M_solidVel,
                           UInt(0),
                           M_solver->FSIOper()->dFESpace().dof().numTotalDof() );
}

void
MS_Model_FSI3D::updateBCSegregated(const std::string& dataFileName)
{

    M_solver->setup(/*data_file*/);
    /*update bc with bcInterface*/
    M_fluidBC.reset(             new BCInterface<FSIOperator>() );
    M_harmonicExtensionBC.reset( new BCInterface<FSIOperator>() );
    M_linearizedFluidBC.reset(   new BCInterface<FSIOperator>() );
    M_solidBC.reset(             new BCInterface<FSIOperator>() );
    M_linearizedSolidBC.reset(   new BCInterface<FSIOperator>() );

    M_fluidBC->SetOperator(             M_solver->FSIOper() );
    M_harmonicExtensionBC->SetOperator( M_solver->FSIOper() );
    M_linearizedFluidBC->SetOperator(   M_solver->FSIOper() );
    M_solidBC->SetOperator(             M_solver->FSIOper() );
    M_linearizedSolidBC->SetOperator(   M_solver->FSIOper() );
    if ( M_solver->isFluid() )
    {
        M_fluidBC->FillHandler( dataFileName, "fluid" );
        M_harmonicExtensionBC->FillHandler( dataFileName, "mesh_motion" );
        M_linearizedFluidBC->FillHandler( dataFileName, "lin_fluid" );
    }
    if ( M_solver->isSolid() )
    {
        M_solidBC->FillHandler( dataFileName, "solid" );
        M_linearizedSolidBC->FillHandler( dataFileName, "lin_solid" );
    }
    M_solver->setFluidBC(                         M_fluidBC->GetHandler() );
    M_solver->setHarmonicExtensionBC( M_harmonicExtensionBC->GetHandler() );
    M_solver->setLinFluidBC(            M_linearizedFluidBC->GetHandler() );
    M_solver->setSolidBC(                         M_solidBC->GetHandler() );
    M_solver->setLinSolidBC(            M_linearizedSolidBC->GetHandler() );

}

void
MS_Model_FSI3D::updateBCMonolithic(const std::string& dataFileName)
{
    /*update bc without bcInterface*/
//     M_solver->setFluxBC(BCh_monolithicFlux());
//     M_solver->setup(/*data_file*/);
//     M_solver->setRobinBC(BCh_monolithicRobin(*M_solver->FSIOper()));
//     M_solver->setFluidBC(BCh_monolithicFluid(*M_solver->FSIOper()));
//     M_solver->setHarmonicExtensionBC (BCh_harmonicExtension(*M_solver->FSIOper()));
//     M_solver->setSolidBC(BCh_monolithicSolid(*M_solver->FSIOper()));

    /*update bc with bcInterface (preliminary version)*/
    M_solver->setup(/*data_file*/);

    M_fluidBC.reset(             new BCInterface<FSIOperator>() );
    M_harmonicExtensionBC.reset( new BCInterface<FSIOperator>() );
    M_solidBC.reset(             new BCInterface<FSIOperator>() );
    M_fluidBC->SetOperator(             M_solver->FSIOper() );
    M_harmonicExtensionBC->SetOperator( M_solver->FSIOper() );
    //M_solidBC->SetOperator(             M_solver->FSIOper() );
    M_fluidBC->FillHandler( dataFileName, "fluid" );
    M_harmonicExtensionBC->FillHandler( dataFileName, "mesh_motion" );
    M_solidBC->FillHandler( dataFileName, "solid" );
    M_solver->setFluidBC(
                   M_fluidBC->GetHandler() );
    M_solver->setHarmonicExtensionBC( M_harmonicExtensionBC->GetHandler() );
    M_solver->setSolidBC(                         M_solidBC->GetHandler() );
}

} // Namespace LifeV
