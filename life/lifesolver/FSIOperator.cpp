/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
#ifndef TWODIM
#include <life/lifesolver/FSIOperator.hpp>
//#include <life/lifesolver/reducedLinFluid.hpp>

namespace LifeV {

// ===================================================
//! Constructors & Destructors
// ===================================================
FSIOperator::FSIOperator():
    M_uFESpace                           ( ),
    M_pFESpace                           ( ),
    M_dFESpace                           ( ),
    M_mmFESpace                          ( ),
    M_fluidMeshPart                      ( ),
    M_solidMeshPart                      ( ),
    M_BCh_u                              ( new fluid_bchandler_raw_type ),
    M_BCh_d                              ( new solid_bchandler_raw_type ),
    M_BCh_mesh                           ( new fluid_bchandler_raw_type ),
    M_BCh_du                             ( new fluid_bchandler_raw_type ),
    M_BCh_du_inv                         ( new fluid_bchandler_raw_type ),
    M_BCh_dz                             ( new solid_bchandler_raw_type ),
    M_BCh_dz_inv                         ( new solid_bchandler_raw_type ),
    M_BCh_dp                             ( new BCHandler ),
    M_BCh_dp_inv                         ( new BCHandler ),
    M_fluid                              ( ),
    M_solid                              ( ),
    M_meshMotion                         ( ),
//     M_fluidLin                           ( ),
//     M_solidLin                           ( ),
    M_bdf                                ( ),
    M_dataFile                           ( ),
    M_dataFluid                          ( new data_fluid() ),
    M_dataSolid                          ( ),
    M_fluidInterfaceMap                  ( ),
    M_solidInterfaceMap                  ( ),
    M_fluidInterfaceMapOnZero            ( ),
    M_solidInterfaceMapOnZero            ( ),
//     M_dofInterfaceFluid                  ( ),
//     M_dofInterfaceSolid                  ( ),
    M_dofFluidToStructure                ( new DofInterface3Dto3D ),
//     M_dofFluidToSolid                    ( new DofInterface3Dto3D ),
//     M_dofSolidToFluid                    ( new DofInterface3Dto3D ),
    M_dofStructureToFluid                ( new DofInterface3Dto3D ),
    M_dofStructureToSolid                ( new DofInterface3Dto3D ),
    M_dofStructureToHarmonicExtension    ( new DofInterface3Dto3D ),
    M_dofHarmonicExtensionToFluid        ( new DofInterface3Dto3D ),
//     M_dofStructureToReducedFluid         ( new DofInterface3Dto3D ),
//     M_dofReducedFluidToStructure         ( new DofInterface3Dto3D ),
    M_fluidInterfaceFlag                 ( 1 ),
    M_solidInterfaceFlag                 ( 1 ),
    M_structureInterfaceFlag             ( 1 ),
    M_harmonicInterfaceFlag              ( 1 ),
    M_interfaceTolerance                 ( 0. ),
    M_dofFluid                           ( new DofInterface3Dto2D ),
    M_dofSolid                           ( new DofInterface3Dto2D ),
    M_dofFluidInv                        ( new DofInterface3Dto2D ),
    M_dofSolidInv                        ( new DofInterface3Dto2D ),
    M_bcvFluidInterfaceDisp              ( new  BCVectorInterface ),
    M_bcvFluidLoadToStructure            ( new  BCVectorInterface ),
    M_bcvSolidLoadToStructure            ( new  BCVectorInterface ),
    M_bcvStructureToFluid                ( new  BCVectorInterface ),
    M_bcvStructureDispToFluid            ( new  BCVectorInterface ),
    M_bcvStructureDispToSolid            ( new  BCVectorInterface ),
    M_bcvStructureDispToHarmonicExtension( new  BCVectorInterface ),
    M_bcvHarmonicExtensionVelToFluid     ( new  BCVectorInterface ),
//     M_bcvStructureToReducedFluid         ( new  BCVectorInterface ),
//     M_bcvReducedFluidToStructure         ( new  BCVectorInterface ),
    M_bcvDerHarmonicExtensionVelToFluid  ( new  BCVectorInterface ),
    M_bcvDerFluidLoadToStructure         ( new  BCVectorInterface ),
    M_bcvDerFluidLoadToFluid             ( new  BCVectorInterface ),
    M_bcvDerStructureDispToSolid         ( new  BCVectorInterface ),
//     M_bcvDerReducedFluidLoadToStructure  ( new  BCVectorInterface ),
//     M_bcvDerStructureAccToReducedFluid   ( new  BCVectorInterface ),
    M_lambdaFluid                        ( ),
    M_lambdaFluidRepeated                ( ),
    M_lambda                             ( ),
    M_lambdaDot                             ( ),
    M_lambdaSolid                        ( ),
    M_lambdaSolidRepeated                ( ),
    M_lambdaSolidOld                     ( ),
    M_lambdaDotSolid                     ( ),
    M_lambdaDotSolidRepeated             ( ),
    M_sigmaFluid                         ( ),
    M_sigmaSolid                         ( ),
    M_sigmaFluidRepeated                 ( ),
    M_sigmaSolidRepeated                 ( ),
    M_minusSigmaFluid                    ( ), //sigmafluid: auxiliary variable for Robin.
    M_minusSigmaFluidRepeated            ( ), //sigmafluid: auxiliary variable for Robin.
    M_dispFluidMeshOld                   ( ),
    M_veloFluidMesh                      ( ),
    M_derVeloFluidMesh                   ( ),
    M_un                                 ( ),
    M_rhs                                ( ),
    M_Alphaf                             ( ), //vector_type, for alphaf robin
    M_AlphafCoef                         ( 0 ),
    M_betamedio                          ( ),
    M_nbEval                             ( 0 ),
    M_epetraComm                         ( ),
    M_epetraWorldComm                    ( ),
    M_method                             ( ),
    M_algorithm                          ( ),
    M_precond                            ( NO_PRECONDITIONER ),
    M_DDNprecond                         ( ),
    M_mpi                                ( true ),
    M_isFluid                            ( false ),
    M_isSolid                            ( false ),
    M_linearFluid                        ( ),
    M_linearSolid                        ( ),
    M_fluidLeader                        ( ),
    M_solidLeader                        ( )
{
}

FSIOperator::~FSIOperator()
{
}




// ===================================================
//! Virtual Methods
// ===================================================
void
FSIOperator::setDataFromGetPot( const GetPot& dataFile )
{
	M_dataFile               = dataFile;

    M_dataFluid->setup( M_dataFile );
    M_dataSolid.reset( new data_solid(dataFile) );

    M_method                 = dataFile( "problem/method",    "steklovPoincare" );
    M_algorithm              = dataFile( "problem/algorithm", "DirichletNeumann" );

    M_fluidInterfaceFlag     = dataFile( "interface/fluid_flag",     M_fluidInterfaceFlag );
    M_solidInterfaceFlag     = dataFile( "interface/solid_flag",     M_fluidInterfaceFlag );
    M_structureInterfaceFlag = dataFile( "interface/structure_flag", M_fluidInterfaceFlag );
    M_harmonicInterfaceFlag  = dataFile( "interface/harmonic_flag",  M_fluidInterfaceFlag );
    M_interfaceTolerance     = dataFile( "interface/tolerance",      0. );

	M_epetraWorldComm->Barrier();
}



void
FSIOperator::setupFEspace()
{
    if( M_epetraComm->MyPID()==0)
        std::cout<< "FSIOperator: setting RefFE and QuadRule ... \n";
	std::string uOrder = M_dataFluid->uOrder();
	std::string pOrder = M_dataFluid->pOrder();
	std::string dOrder = M_dataSolid->order();

    const RefFE*    refFE_vel(0);
    const QuadRule* qR_vel(0);
    const QuadRule* bdQr_vel(0);

    const RefFE*    refFE_press(0);
    const QuadRule* qR_press(0);
    const QuadRule* bdQr_press(0);

    const RefFE*    refFE_struct(0);
    const QuadRule* qR_struct(0);
    const QuadRule* bdQr_struct(0);

    Displayer disp(M_epetraComm.get());
    disp.leaderPrint("velocity order = ", uOrder, "\n");
    if ( uOrder.compare("P2") == 0 )
    {
        refFE_vel = &feTetraP2;
        qR_vel    = &quadRuleTetra15pt; // DoE 5
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else
    	if ( uOrder.compare("P1") == 0 )
    	{
			refFE_vel = &feTetraP1;
            qR_vel    = &quadRuleTetra4pt;  // DoE 2
            bdQr_vel  = &quadRuleTria3pt;   // DoE 2
        }
        else
            if ( uOrder.compare("P1Bubble") == 0 )
            {
                refFE_vel = &feTetraP1bubble;
                qR_vel    = &quadRuleTetra64pt;  // DoE 2
                bdQr_vel  = &quadRuleTria3pt;   // DoE 2
            }
            else
            {
                ERROR_MSG(uOrder + " velocity FE not implemented yet.");
            }

    disp.leaderPrint("pressure order = ", pOrder, "\n");
    if ( pOrder.compare("P2") == 0 )
    {
        refFE_press = &feTetraP2;
        qR_press    = qR_vel; // DoE 5
        bdQr_press  = &quadRuleTria3pt;   // DoE 2
    }
    else
    	if ( pOrder.compare("P1") == 0 )
		{
    		refFE_press = &feTetraP1;
//             qR_press    = &quadRuleTetra64pt;  // DoE 2
    		qR_press    = qR_vel;  // DoE 2
    		bdQr_press  = &quadRuleTria3pt;   // DoE 2
    	}
    	else
    	{
    		ERROR_MSG(pOrder +" pressure FE not implemented yet.");
    	}

    disp.leaderPrint("structure order = ", dOrder, "\n");
    if ( dOrder.compare("P2") == 0 )
    {
        refFE_struct = &feTetraP2;
        qR_struct    = &quadRuleTetra15pt; // DoE 5
        bdQr_struct  = &quadRuleTria3pt;   // DoE 2
    }
    else
    	if ( dOrder.compare("P1") == 0 )
    	{
    		refFE_struct = &feTetraP1;
    		qR_struct    = &quadRuleTetra4pt;  // DoE 2
    		bdQr_struct  = &quadRuleTria3pt;   // DoE 2
    	}
    	else
		{
			ERROR_MSG(dOrder + " structure FE not implemented yet.");
		}


    disp.leaderPrint("FSIOperator: building the fluid FESpace ... ");
    if (this->isFluid())
    {
        M_fluidMeshPart.reset(new  partitionMesh< mesh_type > (*M_dataFluid->dataMesh()->mesh(), *M_epetraComm));

        M_mmFESpace.reset(new FESpace<mesh_type, EpetraMap>(*M_fluidMeshPart,
															//dOrder,
															*refFE_struct,
                                                            *qR_struct,
                                                            *bdQr_struct,
                                                            3,
                                                            *M_epetraComm));

        M_uFESpace.reset( new FESpace<mesh_type, EpetraMap>(*M_fluidMeshPart,
															//uOrder,
                                                            *refFE_vel,
                                                            *qR_vel,
                                                            *bdQr_vel,
                                                            3,
                                                            *M_epetraComm));

        M_pFESpace.reset( new FESpace<mesh_type, EpetraMap>(*M_fluidMeshPart,
															//pOrder,
                                                            *refFE_press,
                                                            *qR_press,
                                                            *bdQr_press,
                                                            1,
                                                            *M_epetraComm));
    }
    else
    {
        M_mmFESpace.reset(new FESpace<mesh_type, EpetraMap>(M_dataFluid->dataMesh()->mesh(),
															//dOrder,
															*refFE_struct,
                                                            *qR_struct,
                                                            *bdQr_struct,
                                                            3,
                                                            *M_epetraComm));

        M_uFESpace.reset( new FESpace<mesh_type, EpetraMap>(M_dataFluid->dataMesh()->mesh(),
															//uOrder,
                                                            *refFE_vel,
                                                            *qR_vel,
                                                            *bdQr_vel,
                                                            3,
                                                            *M_epetraComm));

        M_pFESpace.reset( new FESpace<mesh_type, EpetraMap>(M_dataFluid->dataMesh()->mesh(),
															//pOrder,
                                                            *refFE_press,
                                                            *qR_press,
                                                            *bdQr_press,
                                                            1,
                                                            *M_epetraComm));
    }
    M_epetraWorldComm->Barrier();
    disp.leaderPrint("fluid: ok.\n");





    disp.leaderPrint("FSIOperator: building the solid FESpace ... ");
    if (this->isSolid())
    {
    	M_solidMeshPart.reset( new  partitionMesh< mesh_type > ( *M_dataSolid->mesh(), *M_epetraComm ) );
    	M_dFESpace.reset( new FESpace<mesh_type, EpetraMap>( *M_solidMeshPart,
    															 dOrder,
    	                                                         //*refFE_struct,
    	                                                         //*qR_struct,
    	                                                         //*bdQr_struct,
    	                                                         3,
    	                                                         *M_epetraComm));
    }
    else
    {
        M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(M_dataSolid->mesh(),
                                                           //dOrder,
                                                           *refFE_struct,
                                                           *qR_struct,
                                                           *bdQr_struct,
                                                           3,
                                                           *M_epetraComm));
    }
    M_epetraWorldComm->Barrier();
	disp.leaderPrint("solid: ok.\n");
}



void
FSIOperator::setupDOF( void )
{
    Displayer disp(M_epetraWorldComm.get());
    disp.leaderPrint("FSIOperator: setting DOF ... " );
    Dof uDof(*M_dataFluid->dataMesh()->mesh(), M_uFESpace->refFE());
//     Dof pDof(*M_dataFluid->mesh(), M_pFESpace->refFE());
    Dof dDof(*M_dataSolid->mesh(), M_dFESpace->refFE());

	M_dofFluidToStructure->setup(   M_dFESpace->refFE(), dDof, //M_dFESpace->dof(),
			                        M_uFESpace->refFE(), M_uFESpace->dof() );
	M_dofFluidToStructure->update( *M_dataSolid->mesh(), M_fluidInterfaceFlag,
								   *M_uFESpace->mesh(),  M_structureInterfaceFlag,
								    M_interfaceTolerance );

	//here the solid mesh must be non partitioned in the monolithic case
	M_dofStructureToHarmonicExtension->setup(   M_uFESpace->refFE(), M_uFESpace->dof(),
											    M_dFESpace->refFE(), M_dFESpace->dof() );
	M_dofStructureToHarmonicExtension->update( *M_uFESpace->mesh(),  M_structureInterfaceFlag,
											   *M_dFESpace->mesh(),  M_harmonicInterfaceFlag,
											    M_interfaceTolerance );

	M_dofStructureToSolid->setup(   M_dFESpace->refFE(), M_dFESpace->dof(),
								    M_dFESpace->refFE(), M_dFESpace->dof() );
	M_dofStructureToSolid->update( *M_dFESpace->mesh(),  M_structureInterfaceFlag,
								   *M_dFESpace->mesh(),  M_solidInterfaceFlag,
								    M_interfaceTolerance );

	M_dofStructureToFluid->setup(   M_uFESpace->refFE(), M_uFESpace->dof(), //modifica matteo FSI
								    M_dFESpace->refFE(), M_dFESpace->dof() );
	M_dofStructureToFluid->update( *M_uFESpace->mesh(),  M_structureInterfaceFlag,
								   //*M_dataFluid->mesh(), M_structureInterfaceFlag,
								   *M_dataSolid->mesh(), M_fluidInterfaceFlag,
								    M_interfaceTolerance );

	M_dofHarmonicExtensionToFluid->setup(   M_uFESpace->refFE(),  uDof, //M_uFESpace->dof(),
										    M_uFESpace->refFE(),  uDof); //M_uFESpace->dof() );
	M_dofHarmonicExtensionToFluid->update( *M_dataFluid->dataMesh()->mesh(),  M_harmonicInterfaceFlag,
										   *M_dataFluid->dataMesh()->mesh(),  M_fluidInterfaceFlag,
										    M_interfaceTolerance );

	M_epetraWorldComm->Barrier();
	disp.leaderPrint(" done.\n");

    createInterfaceMaps(M_dofStructureToHarmonicExtension);
}



void FSIOperator::createInterfaceMaps(dof_interface_type3D dofStructureToHarmonicExtension)
{

    Displayer disp(M_epetraWorldComm.get());
	// now we build the sigma and lambda variables on each proc
	disp.leaderPrint("FSIOperator: building fluid variables ... ");

	std::map<ID, ID> const& locDofMap = dofStructureToHarmonicExtension->locDofMap();

	std::vector<int> dofInterfaceFluid;
	dofInterfaceFluid.reserve( M_dofHarmonicExtensionToFluid->locDofMap().size() );

	//is the interface map between HE (first) and solid (second)
	if( this->isFluid() )
	{
		for (int dim = 0; dim < (int)nDimensions; ++dim)
			for ( Iterator i = locDofMap.begin(); i != locDofMap.end(); ++i )
				dofInterfaceFluid.push_back(i->second + dim * M_dFESpace->dof().numTotalDof()); // in solid numerotation
	}

	int* pointerToDofs(0);
	if (dofInterfaceFluid.size() > 0)
		pointerToDofs = &dofInterfaceFluid[0];

	M_fluidInterfaceMap.reset( new EpetraMap( -1,
											  static_cast<int>(dofInterfaceFluid.size()),
											  pointerToDofs,
											  1,
											  *M_epetraWorldComm ));
	disp.leaderPrint(" done.\n");
	M_epetraWorldComm->Barrier();



	disp.leaderPrint("FSIOperator: building solid variables ... ");

	std::vector<int> dofInterfaceSolid;
	dofInterfaceSolid.reserve(M_dofStructureToSolid->locDofMap().size());

	if (this->isSolid())
	{
		//std::cout << "solid" << std::endl;
		for (int dim = 0; dim < (int)nDimensions; ++dim)
			for ( Iterator i = locDofMap.begin(); i != locDofMap.end(); ++i )
					dofInterfaceSolid.push_back(i->second + dim * M_dFESpace->dof().numTotalDof()); // in solid numerotation
	}


	pointerToDofs = 0;
	if (dofInterfaceSolid.size() > 0)
		pointerToDofs = &dofInterfaceSolid[0];

	M_solidInterfaceMap.reset( new EpetraMap( -1,
											  static_cast<int>(dofInterfaceSolid.size()),
										      pointerToDofs,
											  1,
											  *M_epetraWorldComm ));

	M_epetraWorldComm->Barrier();
	disp.leaderPrint(" done.\n");

	//    M_dofStructureToHarmonicExtension->showMe(true, std::cout);
	//    M_dofHarmonicExtensionToFluid->showMe(true, std::cout);
	//    M_dofStructureToReducedFluid->setup(uFESpace.refFE(), M_fluid->pressFESpace().dof(),
	//                                             dFESpace.refFE(), dFESpace.dof());
	//         M_dofStructureToReducedFluid->update(fluidMesh, 1,
	//                                              solidMesh, 1,
	//                                              0.0);

	//         M_dofReducedFluidToStructure->setup(dFESpace.refFE(), dFESpace.dof(),
	//                                             uFESpace.refFE(), uFESpace.dof());
	//         M_dofReducedFluidToStructure->update(solidMesh, 1,
	//                                              fluidMesh, 1,
	//                                              0.0);
	//     }



	disp.leaderPrint("FSIOperator: variables initialization ... ");

	//variablesInit( refFE_struct, bdQr_struct, qR_struct);
	variablesInit( M_dataSolid->order() );

	M_epetraWorldComm->Barrier();
	disp.leaderPrint(" done.\n");

}




void
FSIOperator::setupFluidSolid( void )
{
    if ( this->isFluid() )
    {
        UInt numLM = imposeFlux();

        M_meshMotion.reset( new meshmotion_raw_type(               *M_mmFESpace,             *M_epetraComm ) );
        M_fluid.reset(      new fluid_raw_type(      *M_dataFluid, *M_uFESpace, *M_pFESpace, *M_epetraComm, numLM ) );
        M_solid.reset(      new solid_raw_type(      *M_dataSolid, *M_dFESpace,              *M_epetraComm ) );

//         if ( M_linearFluid )
//             M_fluidLin.reset( new FSIOperator::fluidlin_raw_type( *M_dataFluid, *M_uFESpace, *M_pFESpace, *M_epetraComm ) );

//         if ( M_linearSolid )
//             M_solidLin.reset( new FSIOperator::solidlin_raw_type( *M_dataSolid, *M_dFESpace, *M_epetraComm ) );

        //Vector initialization
    	M_rhs.reset( new vector_type( M_fluid->getMap() ) );
    }

    if ( this->isSolid() )
    {
//         M_fluid.reset( new fluid_raw_type( *M_dataFluid, *M_uFESpace, *M_pFESpace, *M_epetraComm ) );
    	M_meshMotion.reset( new meshmotion_raw_type(               *M_mmFESpace, *M_epetraComm ) );
        M_solid.reset(      new solid_raw_type(      *M_dataSolid, *M_dFESpace,  *M_epetraComm ) );

//         if ( M_linearFluid )
//             M_fluidLin.reset( new FSIOperator::fluidlin_raw_type( *M_dataFluid, *M_uFESpace, *M_pFESpace, *M_epetraComm ) );

//         if ( M_linearSolid )
//             M_solidLin.reset( new FSIOperator::solidlin_raw_type( *M_dataSolid, *M_dFESpace, *M_epetraComm ) );
    }

	M_epetraWorldComm->Barrier();
}



void
FSIOperator::setupSystem( void )
{
    if ( this->isFluid() )
    {
        //Data
        M_fluid->setUp( M_dataFile );
        M_meshMotion->setUp( M_dataFile );
//         if (M_linearFluid)
//             M_fluidLin->setUp(dataFile);
    }

    if ( this->isSolid() )
    {
        M_solid->setUp( M_dataFile );
//         if (M_linearSolid)
//             M_solidLin->setUp( dataFile );
    }

	M_epetraWorldComm->Barrier();
}



void
FSIOperator::buildSystem()
{
    if ( this->isFluid() )
    {
        M_fluid->buildSystem();
//         if (M_linearFluid)
//             M_fluidLin->buildSystem();
    }

    if ( this->isSolid() )
    {
        M_solid->buildSystem();
//         if (M_linearSolid)
//             M_solidLin->buildSystem();
    }

    M_epetraWorldComm->Barrier();
}



void
FSIOperator::updateSystem( )
{
    shiftSolution();

    if ( this->isFluid() )
    {
        M_meshMotion->updateSystem();

        transferMeshMotionOnFluid(M_meshMotion->disp(), *this->M_dispFluidMeshOld);

        if(M_fluid->solution().get())
            M_un                = M_fluid->solution();
        *M_rhs               = M_fluid->matrMass()*M_bdf->time_der( M_dataFluid->dataTime()->getTimeStep() );
    }

    if ( this->isSolid() )
    {
        this->M_solid->updateSystem();
    }
    couplingVariableExtrap();
}



void FSIOperator::couplingVariableExtrap( )
{
	*M_lambda      = lambdaSolid();
	if (!M_lambdaDot.get())
    {
        M_lambdaDot.reset        ( new vector_type(*M_fluidInterfaceMap, Unique) );
        *M_lambda     += M_dataFluid->dataTime()->getTimeStep()*lambdaDotSolid();
    }
    else
    {
        *M_lambda     += 1.5*M_dataFluid->dataTime()->getTimeStep()*lambdaDotSolid(); // *1.5
        *M_lambda     -= M_dataFluid->dataTime()->getTimeStep()*0.5*(*M_lambdaDot);
    }

	*M_lambdaDot   = lambdaDotSolid();

	displayer().leaderPrint("\n norm( disp  ) init = ", M_lambda->NormInf() );
	displayer().leaderPrint("\n norm( velo )  init = ", M_lambdaDot->NormInf());
}



void
FSIOperator::shiftSolution()
{
    if ( this->isFluid() )
    {
        this->M_bdf->shift_right( *M_fluid->solution() );
    }
}



void
FSIOperator::initializeFluid( const vector_type& velAndPressure,
                              const vector_type& displacement )
{
    this->fluid().initialize( velAndPressure );
    this->meshMotion().initialize( displacement );
    this->moveMesh( displacement);
}



void
FSIOperator::initializeSolid( vector_ptrtype displacement,
                              vector_ptrtype velocity )
{
    this->solid().initialize( displacement, velocity);
}



void
FSIOperator::updateJacobian( const vector_type& /*sol*/, const int& /*iter*/)
{
}



void
FSIOperator::moveMesh( const vector_type& dep )
{
	displayer().leaderPrint( "  Moving the mesh ... ");
    M_fluidMeshPart->mesh()->moveMesh(dep,  this->M_mmFESpace->dof().numTotalDof());
    displayer().leaderPrint(  " done.\n" );
    M_fluid->recomputeMatrix(true);
}



void
FSIOperator::transferFluidOnInterface(const vector_type &_vec1, vector_type &_vec2)
{
    // e.g.: vec1=M_fluid->residual(), vec2=M_sigmaFluid
//     _vec2 = ZeroVector(_vec2.size());

    if (_vec1.getMaptype() == Unique)
	{
		vector_type const  vec1Repeated(_vec1, Repeated);
		transferFluidOnInterface(vec1Repeated, _vec2);
		return;
	}

    if (_vec2.getMaptype() == Repeated)
	{
		vector_type  vec2Unique(_vec2, Unique);
		transferFluidOnInterface(_vec1, vec2Unique);
		_vec2 = vec2Unique;
		return;
	}

    _vec2 *= 0;

    std::map<ID, ID> const& locDofMap = M_dofFluidToStructure->locDofMap();


    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();
    int numTotalDofFluid = M_uFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::const_iterator Iterator;

    for (int dim = 0; dim < (int)nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
		{
//                 std::cout <<  it->second + dim*numTotalDofFluid << " to "
//                           <<  it->first + dim*numTotalDofSolid << " : "
//                           << _vec1[it->second + dim*numTotalDofFluid] << std::endl;
			_vec2.checkAndSet( it->first + dim*numTotalDofSolid,
							   _vec1[it->second + dim*numTotalDofFluid] );
		}
}



//works in serial but no yet in parallel
void
FSIOperator::transferSolidOnFluid(const vector_type &_vec1, vector_type &_vec2)
{
    //    e.g.: vec2=M_fluid->residual(), vec1=M_sigmaFluid
    //     _vec2 = ZeroVector(_vec2.size());

    if (_vec1.getMaptype() == Unique)
	{
		vector_type const  vec1Repeated(_vec1, Repeated);
		transferSolidOnInterface(vec1Repeated, _vec2);
		return;
	}

    if (_vec2.getMaptype() == Repeated)
	{
		vector_type  vec2Unique(_vec2, Unique);
		transferSolidOnInterface(_vec1, vec2Unique);
		_vec2 = vec2Unique;
		return;
	}

    _vec2 *= 0;

    std::map<ID, ID> const& locDofMap = M_dofFluidToStructure->locDofMap();


    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();
    int numTotalDofFluid = M_uFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::const_iterator Iterator;

    for (UInt dim = 0; dim < nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
		{
			_vec2.checkAndSet( it->second + dim*numTotalDofFluid,
							   _vec1[it->first + dim*numTotalDofSolid]);
		}
}



void
FSIOperator::transferSolidOnInterface(const vector_type &_vec1, vector_type &_vec2)
{
    /* e.g.:
       vec1 (Unique)       vec2 (Repeated)  (On different EpetraMaps) (Changing now to Unique on both)
       M_solid->disp()     M_lambdaSolid
       M_solid->vel()      M_lambdaDotSolid
       M_solid->residual() M_sigmaSolid
    */

    if (_vec1.getMaptype() == Unique)
	{
		vector_type const  vec1Repeated(_vec1, Repeated);
		transferSolidOnInterface(vec1Repeated, _vec2);
		return;
	}

    if (_vec2.getMaptype() == Repeated)
	{
		vector_type  vec2Unique(_vec2, Unique);
		transferSolidOnInterface(_vec1, vec2Unique);
		_vec2 = vec2Unique;
		return;
	}

    _vec2 *= 0;

    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();

    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::const_iterator Iterator;

    for (UInt dim = 0; dim < nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
		{
			_vec2.checkAndSet( it->second + dim*numTotalDofSolid,
							   _vec1[it->second + dim*numTotalDofSolid] );
		}
}

void
FSIOperator::transferInterfaceOnSolid(const vector_type& _vec1, vector_type& _vec2)
{
    /* e.g.:
       vec2                vec1
       M_solid->disp()     M_lambdaSolid
       M_solid->vel()      M_lambdaDotSolid
       M_solid->residual() M_sigmaSolid
    */

    if (_vec1.getMaptype() == Unique)
	{
		vector_type const  vec1Repeated(_vec1, Repeated);
		transferInterfaceOnSolid(vec1Repeated, _vec2);
		return;
	}

    if (_vec2.getMaptype() == Repeated)
	{
		vector_type  vec2Unique(_vec2, Unique);
		transferInterfaceOnSolid(_vec1, vec2Unique);
		_vec2 = vec2Unique;
		return;
	}

    _vec2 *= 0;

    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::const_iterator Iterator;

    for (int dim = 0; dim < (int)nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
		{
			_vec2.checkAndSet( it->second + dim*numTotalDofSolid,
							   _vec1[it->second + dim*numTotalDofSolid] );
		}
}





// ===================================================
//! Display Methods
// ===================================================
bool
FSIOperator::isLeader() const
{
    if ( isFluid() )
    {
        if ( M_fluid.get() == 0 )
            return ( M_epetraComm->MyPID() == 0 );

        return M_fluid->getDisplayer().isLeader();
    }

    if ( M_solid.get() == 0 )
        return ( M_epetraComm->MyPID() == 0 );

    return M_solid->getDisplayer().isLeader();
}



Displayer const&
FSIOperator::displayer()
{
    if ( isFluid() &&  M_fluid.get())
        return M_fluid->getDisplayer();

    if( !isSolid() || !M_solid.get() )
        std::cout << "ERROR: displayer not ready" << std::endl;

    return M_solid->getDisplayer();
}





// ===================================================
//! Get Functions
// ===================================================
/*
FSIOperator::vector_type
FSIOperator::displacementOnInterface()
{

//     vector_type dispOnInterface();
// //@     dispOnInterface = ZeroVector(dispOnInterface.size());

//  //    FOR_EACH_INTERFACE_DOF( dispOnInterface[IDsolid - 1 + jDim*totalDofSolid] =
// //                             M_solid->disp()[IDsolid - 1 + jDim*totalDofSolid]);

//     Real norminf;
//     Real norm2;

//     dispOnInterface.NormInf(&norminf);
//     dispOnInterface.Norm2(&norm2);

//     std::cout << "max norm disp = " << norminf;
//     std::cout << std::endl;
//     std::cout << "l2  norm disp = " << norm2;
//     std::cout << std::endl;

//     return dispOnInterface;
    assert(false);

}
*/





// ===================================================
//! Set Functions
// ===================================================
void
FSIOperator::setComm( const boost::shared_ptr<Epetra_MpiComm>& comm,
                      const boost::shared_ptr<Epetra_MpiComm>& worldComm )
{
    M_epetraComm       = comm;
    M_epetraWorldComm  = worldComm;
}



void
FSIOperator::setTime( const Real& time )
{
	M_time = time;
	M_dataFluid->dataTime()->setTime(time);
	M_dataSolid->setTime(time);
}



void
FSIOperator::setFluid( const fluid_type& fluid, const meshmotion_type& meshmotion )
{
	M_fluid = fluid;
	M_meshMotion = meshmotion;
	M_isFluid = true;
}



void
FSIOperator::setSolid( const solid_type& solid )
{
	M_solid = solid;
	M_isSolid = true;
}



void
FSIOperator::setFluidBC( const fluid_bchandler_type& bc_fluid )
{
    if ( isFluid() )
	{
		M_BCh_u = bc_fluid;
	}
}



void
FSIOperator::setHarmonicExtensionBC( const fluid_bchandler_type& bc_he )
{
    if ( isFluid() )
	{
		M_BCh_mesh = bc_he;
		//M_meshMotion->setBC(*M_BCh_mesh);
	}
}



void
FSIOperator::setSolidBC( const solid_bchandler_type& bc_solid )
{
    if ( isSolid() )
	{
		M_BCh_d = bc_solid;
	}
}



void
FSIOperator::setLambdaFluid( const vector_type& lambda )
{
    if ( lambda.getMaptype() == Unique )
        *M_lambdaFluid = lambda;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_lambdaFluidRepeated = *M_lambdaFluid;
}



void
FSIOperator::setLambdaSolid( const vector_type& lambda )
{
    if ( lambda.getMaptype() == Unique )
        *M_lambdaSolid = lambda;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_lambdaSolidRepeated = *M_lambdaSolid;
}



void
FSIOperator::setLambdaSolidOld( const vector_type& lambda )
{
    if ( lambda.getMaptype() == Unique )
        *M_lambdaSolidOld = lambda;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry
}



void
FSIOperator::setLambdaDotSolid( const vector_type& lambda )
{
    if ( lambda.getMaptype() == Unique )
        *M_lambdaDotSolid = lambda;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_lambdaDotSolidRepeated = *M_lambdaDotSolid;
}



void
FSIOperator::setSigmaFluid( const vector_type& sigma )
{
    if ( sigma.getMaptype() == Unique )
        *M_sigmaFluid = sigma;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_sigmaFluidRepeated = *M_sigmaFluid;

}


void
FSIOperator::setSigmaSolid( const vector_type& sigma )
{
    if ( sigma.getMaptype() == Unique )
        *M_sigmaSolid = sigma;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_sigmaSolidRepeated = *M_sigmaSolid;
}



void
FSIOperator::setMinusSigmaFluid( const vector_type& sigma )
{
    if ( sigma.getMaptype() == Unique )
    {
    	*M_minusSigmaFluid  = sigma;
    	*M_minusSigmaFluid *= -1;
    }
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_minusSigmaFluidRepeated = *M_minusSigmaFluid;
}



void
FSIOperator::setAlphafCoef( )
{
    Real h=0.1, R=0.5;

    M_AlphafCoef  = 2*(this->dataSolid().rho()*h)/this->dataFluid().dataTime()->getTimeStep();
    M_AlphafCoef += h*this->dataSolid().young(0)*this->dataFluid().dataTime()->getTimeStep() /
                    (2*pow(R,2) *(1-pow(dataSolid().poisson(0),2)));
}



void
FSIOperator::setAlphafbcf( const bc_function_type& alphafbcf )
{
	vector_type vec( M_fluid->velFESpace().map());
	M_fluid->velFESpace().interpolate(alphafbcf, vec, 0.0);
	*M_Alphaf = vec ;
}



void
FSIOperator::setStructureToFluidParametres()
{
	this->setAlphafCoef();
	this->setAlphaf();

	if(M_Alphaf.get()==0)
	{
		this->setAlphafCoef();
		M_bcvStructureToFluid->setMixteCoef(M_AlphafCoef);
		M_bcvStructureToFluid->setBetaCoef(M_AlphafCoef);
	}
	else
	{
		M_bcvStructureToFluid->setMixteVec(this->Alphaf());
		M_bcvStructureToFluid->setBetaVec(this->Alphaf());
	}
}



void
FSIOperator::setStructureDispToHarmonicExtension( const vector_type& disp, UInt type )
{
    M_bcvStructureDispToHarmonicExtension->setup( disp,
                                                  M_dFESpace->dof().numTotalDof(),
                                                  M_dofStructureToHarmonicExtension,
                                                  type );
}



void
FSIOperator::setStructureToFluid( const vector_type& velo,  UInt type )
{
    M_bcvStructureToFluid->setup( velo,
                                  M_uFESpace->dof().numTotalDof(),
                                  M_dofHarmonicExtensionToFluid,
                                  type );
}



void
FSIOperator::setStructureDispToFluid( const vector_type& disp,  UInt type )
{
  M_bcvStructureDispToFluid->setup( disp,
                                    M_uFESpace->dof().numTotalDof(),
				                    M_dofStructureToFluid,
				                    type );
}



void
FSIOperator::setStructureDispToSolid( const vector_type& disp, UInt type )
{
    M_bcvStructureDispToSolid->setup( disp,
                                      M_dFESpace->dof().numTotalDof(),
                                      M_dofStructureToSolid,
                                      type );
}



void
FSIOperator::setDerStructureDispToSolid( const vector_type& ddisp, UInt type )
{
    M_bcvDerStructureDispToSolid->setup( ddisp,
                                         M_dFESpace->dof().numTotalDof(),
                                         M_dofStructureToSolid,
                                         type );
}



void
FSIOperator::setSolidLoadToStructure( const vector_type& load, UInt type )
{
    M_bcvSolidLoadToStructure->setup( load,
                                      M_dFESpace->dof().numTotalDof(),
                                      M_dofStructureToFluid,
                                      type );
}



void
FSIOperator::setHarmonicExtensionVelToFluid( const vector_type& vel, UInt type )
{
    M_bcvHarmonicExtensionVelToFluid->setup( vel,
                                             M_uFESpace->dof().numTotalDof(),
                                             M_dofHarmonicExtensionToFluid,
                                             type );
}



void
FSIOperator::setDerHarmonicExtensionVelToFluid( const vector_type& dvel, UInt type )
{
    M_bcvDerHarmonicExtensionVelToFluid->setup( dvel,
                                                M_uFESpace->dof().numTotalDof(),
                                                M_dofHarmonicExtensionToFluid,
                                                type );
}




void
FSIOperator::setFluidLoadToStructure( const vector_type& load, UInt type )
{
    M_bcvFluidLoadToStructure->setup( load,
                                      M_dFESpace->dof().numTotalDof(),
                                      M_dofStructureToSolid,
                                      type );
}



void
FSIOperator::setDerFluidLoadToStructure( const vector_type& dload, UInt type )
{
    M_bcvDerFluidLoadToStructure->setup( dload,
                                         M_dFESpace->dof().numTotalDof(),
                                         M_dofStructureToSolid,
                                         type );
}



void
FSIOperator::setDerFluidLoadToFluid( const vector_type& dload, UInt type )
{
    M_bcvDerFluidLoadToFluid->setup( dload,
                                     M_uFESpace->dof().numTotalDof(),
                                     M_dofHarmonicExtensionToFluid,
                                     type );
}

void FSIOperator::setMixteOuterWall(function_type const& dload, function_type const& E)
{
    M_bcfMixteOuterWall.setFunctions_Mixte(dload,
                                           E);
}


// void
// FSIOperator::setDerReducedFluidLoadToStructure(Vector &dload, UInt type )
// {
//     M_bcvDerReducedFluidLoadToStructure->setup( dload,
//                                                 M_uFESpace->dof().numTotalDof(),
//                                                 M_dofReducedFluidToStructure,
//                                                 type );
// }



// void
// FSIOperator::setDerStructureAccToReducedFluid(Vector &acc, UInt type )
// {
//     M_bcvDerStructureAccToReducedFluid->setup( acc,
//                                                M_dofStructureToReducedFluid,
//                                                type );
// }



// void
// FSIOperator::setReducedLinFluidBC(fluid_bchandler_type bc_dredfluid )
// {
//     M_BCh_dp     = bc_dredfluid;
//     M_reducedLinFluid->setUpBC(bc_dredfluid);
// }



// void
// FSIOperator::setInvReducedLinFluidBC(fluid_bchandler_type bc_invdredfluid )
// {
//     M_BCh_dp_inv = bc_invdredfluid;
// }

//





// ===================================================
//! Protected Methods
// ===================================================
UInt
FSIOperator::imposeFlux( void )
{
    std::vector<BCName> fluxVector = M_BCh_u->getBCWithType( Flux );
    UInt numLM = static_cast<UInt>( fluxVector.size() );
    if( M_epetraComm->MyPID()==0)
        std::cout<< " numLM = "<< numLM<<std::endl;

    UInt offset = M_uFESpace->map().getMap(Unique)->NumGlobalElements()
                + M_pFESpace->map().getMap(Unique)->NumGlobalElements();

    for ( UInt i = 0; i < numLM; ++i )
    	M_BCh_u->setOffset( fluxVector[i], offset + i );

    return numLM;
}



void
FSIOperator::variablesInit( const std::string& /*dOrder*/ )
//FSIOperator::variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct)
{
    M_lambdaFluid.reset        ( new vector_type(*M_fluidInterfaceMap, Unique) );
    M_lambda.reset             ( new vector_type(*M_solidInterfaceMap, Unique) );
    M_lambdaDot.reset             ( new vector_type(*M_solidInterfaceMap, Unique) );
    M_lambdaFluidRepeated.reset( new vector_type(*M_fluidInterfaceMap, Repeated) );

    if ( this->isFluid() )
	{
		M_dispFluidMeshOld.reset( new vector_type(M_uFESpace->map(), Repeated) );
		M_veloFluidMesh.reset   ( new vector_type(M_uFESpace->map(), Repeated) );
		M_Alphaf.reset          ( new vector_type(M_uFESpace->map(), Repeated));

		if ( M_linearFluid )
			M_derVeloFluidMesh.reset( new vector_type(this->M_uFESpace->map(), Repeated) );
	}

    M_sigmaFluid.reset              ( new vector_type(*M_fluidInterfaceMap, Unique) );
    M_sigmaFluidRepeated.reset      ( new vector_type(*M_fluidInterfaceMap, Repeated) );
    M_minusSigmaFluid.reset         ( new vector_type(*M_fluidInterfaceMap, Unique) );
    M_minusSigmaFluidRepeated.reset ( new vector_type(*M_fluidInterfaceMap, Repeated) );

    M_lambdaSolid.reset   ( new vector_type(*M_solidInterfaceMap, Unique) );
    M_lambdaSolidOld.reset( new vector_type(*M_solidInterfaceMap, Unique) );
    M_lambdaDotSolid.reset( new vector_type(*M_solidInterfaceMap, Unique) );
    M_sigmaSolid.reset    ( new vector_type(*M_solidInterfaceMap, Unique) );

    M_lambdaSolidRepeated.reset   ( new vector_type(*M_solidInterfaceMap, Repeated) );
    M_lambdaDotSolidRepeated.reset( new vector_type(*M_solidInterfaceMap, Repeated) );
    M_sigmaSolidRepeated.reset    ( new vector_type(*M_solidInterfaceMap, Repeated) );
}

void FSIOperator::initializeBDF(vector_ptrtype un)
{
    M_un = un;
    //BDF initialization
    M_bdf.reset( new BdfT<vector_type>( M_dataFluid->dataTime()->getBDF_order() ) );
    M_bdf->initialize_unk( *un );
}


void
FSIOperator::transferMeshMotionOnFluid( const vector_type& _vec1, vector_type& _vec2 )
{
    //transferMeshMotionOnFluid should handle the repetition of the interface nodes.
    if (_vec1.getMaptype() == Unique)
	{
		vector_type const  vec1Repeated(_vec1, Repeated);
		transferMeshMotionOnFluid(vec1Repeated, _vec2);
		return;
	}

    if (_vec2.getMaptype() == Repeated)
	{
		vector_type  vec2Unique(_vec2, Unique);
		transferMeshMotionOnFluid(_vec1, vec2Unique);
		_vec2 = vec2Unique;
		return;
	}

    _vec2 *=0;

    interpolateVelocity(_vec1, _vec2);
    return;

    /*

    std::map<ID, ID> const& locDofMap = M_dofFluidToStructure->locDofMap();


    int numTotalDofMesh  = M_mmFESpace->dof().numTotalDof();
    int numTotalDofFluid = M_uFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::iterator Iterator;

    for (int dim = 0; dim < (int)nDimensions; ++dim)
        for ( Iterator it = locDofMap.begin(); it != locDofMap.end(); ++it )
            {
//                 std::cout << " doing: for " << it->second << " to " << it->second
//                           << " _vec1.Map().LID(it->second + dim*numTotalDofMesh) = " << _vec1.Map().LID(it->second + dim*numTotalDofMesh)
//                           << " _vec2.Map().LID(it->second + dim*numTotalDofFluid) = " << _vec2.Map().LID(it->second + dim*numTotalDofFluid)
//                           << std::endl;

                if (_vec1.BlockMap().LID(it->second + dim*numTotalDofMesh) >= 0 )
                {
                    _vec2[it->second + dim*numTotalDofFluid] = _vec1[it->second + dim*numTotalDofMesh];
                }
//                 else
//                 {
//                     std::cout << " not done: from " << it->second << " to " << it->second << std::endl;
//                 }
            }
    */
}

void
FSIOperator::interpolateVelocity( const vector_type& _vec1, vector_type& _vec2 )
{
    assert(_vec1.getMaptype() == Repeated);
    assert(_vec2.getMaptype() == Unique);

    typedef mesh_type::VolumeShape GeoShape; // Element shape

    UInt nDofpV = M_uFESpace->refFE().nbDofPerVertex(); // number of Dof per vertex
    UInt nDofpE = M_uFESpace->refFE().nbDofPerEdge();   // number of Dof per edge
    UInt nDofpF = M_uFESpace->refFE().nbDofPerFace();   // number of Dof per face
    UInt nDofpEl = M_uFESpace->refFE().nbDofPerVolume(); // number of Dof per Volume

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::numFaces;    // Number of element's faces

    //    UInt nDofElem = M_uFESpace->refFE().nbDof; // Number of local dof per element of the M_uFESpace->mesh() (_mesh.getRefFE().nbDof)
    UInt nDofElemMesh = M_mmFESpace->refFE().nbDof();

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's Dof on a Element

    Real x, y, z;
    Vector wLoc( nDofElemMesh * nDimensions );
    ID lDof;

    // Loop on elements of the mesh
    for ( ID iElem = 1; iElem <= M_uFESpace->mesh()->numVolumes(); ++iElem )
    {
        UInt elemId = M_uFESpace->mesh()->volume( iElem ).localId();
        if (elemId != iElem)
            std::cout << " elemId = " << elemId << " iElem = " << iElem << std::endl;

        // Updating the local mesh velocity in this mesh elment
        for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
            for ( ID idof = 0; idof < nDofElemMesh; ++idof )
			{
				wLoc( icmp * nDofElemMesh + idof ) =
					_vec1( icmp * M_mmFESpace->dof().numTotalDof()
						   + M_mmFESpace->dof().localToGlobal( iElem, idof + 1));
			}
        // Vertex based Dof
        if ( nDofpV )
        {

            // loop on element vertices
            for ( ID iVe = 1; iVe <= nElemV; ++iVe )
            {

                // Loop number of Dof per vertex
                for ( ID l = 1; l <= nDofpV; ++l )
                {
                    lDof = ( iVe - 1 ) * nDofpV + l; // Local dof in this element

                    // Nodal coordinates
                    x = M_uFESpace->refFE().xi( lDof - 1 );
                    y = M_uFESpace->refFE().eta( lDof - 1 );
                    z = M_uFESpace->refFE().zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        Real __sum = 0;
                        for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local Dof on the element
                            __sum += wLoc( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi( idof, x, y, z );

                        // Updating interpolated mesh velocity
                        int iDof = icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( iElem, lDof  );
                        _vec2.checkAndSet( iDof ,__sum);

                    }
                }
            }
        }

        // Edge based Dof
        if ( nDofpE )
        {

            // loop on element edges
            for ( ID iEd = 1; iEd <= nElemE; ++iEd )
            {

                // Loop number of Dof per edge
                for ( ID l = 1; l <= nDofpE; ++l )
                {
                    lDof = nDofElemV + ( iEd - 1 ) * nDofpE + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    x = M_uFESpace->refFE().xi( lDof - 1 );
                    y = M_uFESpace->refFE().eta( lDof - 1 );
                    z = M_uFESpace->refFE().zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        Real __sum = 0;
                        for ( ID idof = 0; idof < nDofElemMesh; ++idof )   // Loop on local Dof on the adjacent element
                            __sum += wLoc( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi( idof, x, y, z ); // Problem here with P2

                        // Updating interpolating vector
                        int iDof = icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( iElem, lDof );
                        _vec2.checkAndSet( iDof ,__sum);

                    }
                }
            }
        }

        // Face based Dof
        if ( nDofpF )
        {

            // loop on element faces
            for ( ID iFa = 1; iFa <= nElemF; ++iFa )
            {

                // Loop on number of Dof per face
                for ( ID l = 1; l <= nDofpF; ++l )
                {

                    lDof = nDofElemE + nDofElemV + ( iFa - 1 ) * nDofpF + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    x = M_uFESpace->refFE().xi( lDof - 1 );
                    y = M_uFESpace->refFE().eta( lDof - 1 );
                    z = M_uFESpace->refFE().zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        Real __sum = 0;
                        for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local Dof on the adjacent element
                            __sum += wLoc( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi( idof, x, y, z ); // Problem here with P2

                        // Updating interpolating vector
                        int iDof = icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( iElem, lDof + 1);
                        _vec2.checkAndSet( iDof ,__sum);
                    }
                }
            }
        }

        // Element based Dof
        // Loop on number of Dof per Element
        for ( ID l = 1; l <= nDofpEl; ++l )
        {
            lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

            // Nodal coordinates
            x = M_uFESpace->refFE().xi( lDof - 1 );
            y = M_uFESpace->refFE().eta( lDof - 1 );
            z = M_uFESpace->refFE().zeta( lDof - 1 );

            // Loop on data vector components
            for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
            {

                // Interpolating data at the nodal point
                Real __sum = 0;
                for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local Dof on the adjacent element
                    __sum += wLoc( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi( idof, x, y, z );

                // Updating interpolating vector
//                std::cout << M_uFESpace->dof().localToGlobal( elemId, lDof ) << " ";
//                std::cout << icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( elemId, lDof ) << std::endl;

                int iDof = icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( elemId, lDof );
                _vec2.checkAndSet( iDof, __sum);
            }
        }
    }

}




// this will interpolate dofs values from fespace1 to fespace2
void
FSIOperator::interpolateInterfaceDofs( const FESpace<mesh_type, EpetraMap>& _fespace1,
                                       const vector_type&                   _vec1,
                                       const FESpace<mesh_type, EpetraMap>& _fespace2,
                                       vector_type&                         _vec2,
                                       dof_interface_type3D&                _dofInterface)
{
    assert(_vec1.getMaptype() == Repeated);
    assert(_vec2.getMaptype() == Unique);

    typedef mesh_type::VolumeShape GeoShape; // Element shape


    UInt nDofPerVert1  = _fespace1.refFE().nbDofPerVertex(); // number of Dof per vertex
    UInt nDofPerEdge1  = _fespace1.refFE().nbDofPerEdge();   // number of Dof per edge
    //UInt nDofPerFace1  = _fespace1.refFE().nbDofPerFace;   // number of Dof per face
    //UInt nDofPerElem1  = _fespace1.refFE().nbDofPerVolume; // number of Dof per Volume

    UInt nDofPerVert2  = _fespace2.refFE().nbDofPerVertex(); // number of Dof per vertex
    UInt nDofPerEdge2  = _fespace2.refFE().nbDofPerEdge();   // number of Dof per edge
    //UInt nDofPerFace2  = _fespace2.refFE().nbDofPerFace;   // number of Dof per face
    //UInt nDofPerElem2  = _fespace2.refFE().nbDofPerVolume; // number of Dof per Volume

    //UInt nElemV        = GeoShape::numVertices; // Number of element's vertices
    //UInt nElemE        = GeoShape::numEdges;    // Number of element's edges
    //UInt nElemF        = GeoShape::numFaces;    // Number of element's faces

    //UInt nBFacesVert   = GeoShape::GeoBShape::numVertices;
    //UInt nBFacesEdge   = GeoShape::GeoBShape::numEdges;
    //UInt nBFacesFace   = GeoShape::GeoBShape::numFaces;

    UInt nBEdges1      = _fespace1.mesh()->numBFaces();
    UInt nBEdges2      = _fespace2.mesh()->numBFaces();

    //    UInt nDofElem = M_uFESpace->refFE().nbDof; // Number of local dof per element of the M_uFESpace->mesh() (_mesh.getRefFE().nbDof)
    UInt numTotalDof1  = _fespace1.dof().numTotalDof();
    UInt numTotalDof2  = _fespace2.dof().numTotalDof();

    //UInt nDofElemVert1 = nElemV * nDofPerVert1; // number of vertex's Dof on a Element
    //UInt nDofElemEdge1 = nElemE * nDofPerEdge1; // number of edge's Dof on a Element
    //UInt nDofElemFace1 = nElemF * nDofPerFace1; // number of face's Dof on a Element

    //UInt nDofElemVert2 = nElemV * nDofPerVert2; // number of vertex's Dof on a Element
    //UInt nDofElemEdge2 = nElemE * nDofPerEdge2; // number of edge's Dof on a Element
    //UInt nDofElemFace2 = nElemF * nDofPerFace2; // number of face's Dof on a Element

    UInt numTotalVert1 = _fespace1.mesh()->numGlobalVertices();
    //UInt numTotalEdge1 = _fespace1.mesh()->numGlobalEdges();
    //UInt numTotalFace1 = _fespace1.mesh()->numGlobalFaces();
    //UInt numTotalVol1  = _fespace1.mesh()->numGlobalVolumes();

    UInt numTotalVert2 = _fespace2.mesh()->numGlobalVertices();
    //UInt numTotalEdge2 = _fespace2.mesh()->numGlobalEdges();
    //UInt numTotalFace2 = _fespace2.mesh()->numGlobalFaces();
    //UInt numTotalVol2  = _fespace2.mesh()->numGlobalVolumes();


    std::map<ID, ID> const& locDofMap = _dofInterface->locDofMap();
    std::map<ID, ID>::const_iterator iter;

    //Real x, y, z;
    Real value;


    //    Vector wLoc( nDofElemMesh1 * nDimensions );


    // Loop on elements of the mesh
    if (nDofPerVert1 && nDofPerVert2)
	{
		//            std::cout << "  -> both FESpace have unknowns on their nodes" << std::endl;
		for ( ID iVert = 1; iVert <= _fespace1.mesh()->numVertices(); ++iVert )
		{
			if ((int)_fespace1.mesh()->pointList(iVert).marker() != M_fluidInterfaceFlag) continue;

			ID nodeID = _fespace1.mesh()->pointList(iVert).id();
			// Loop number of Dof per vertex
			for ( ID l = 1; l <= nDofPerVert1; ++l )
			{
				// Loop on data vector components
				for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
				{
					// "Interpolating" data at the nodal point
					iter = locDofMap.find(nodeID);
					Real value = _vec1( icmp*numTotalDof1 + nodeID );
					// now what to what boundary node ( in the solid numerotation ) should we send this value ?
					//std::cout << "" << std::endl;
					int iDof = icmp*numTotalDof2 + iter->second;
					//                                    std::cout << " transfering " << value << " from P1 " << nodeID << " to P1 " << iDof  << " ( " << iter->second << " ) " << std::endl;
					// Updating interpolated mesh velocity
					_vec2.checkAndSet( iDof, value );
				}
			}
		}
	}

    // Edge based Dof
    if ( nDofPerEdge1 )
	{
		// loop on boundary edges
		for ( ID iEdge = 1; iEdge <= nBEdges1; ++iEdge )
		{
			if ((int)_fespace1.mesh()->edgeList(iEdge).marker() != M_fluidInterfaceFlag) continue;

			// edge ID
			ID edgeID = _fespace1.mesh()->edgeList(iEdge).id();
			// dof position of the edge since unknowns are store using [ node | edges | faces | volumes ]
			int iDofEdge = numTotalVert1 + edgeID;

			if (nDofPerEdge2) // the second FE space has dofs on its edges
			{
				for ( ID l = 1; l <= nDofPerEdge1; ++l )
				{
					// Loop on data vector components
					for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
					{
						// ID of the dof in the solid numerotation
						iter = locDofMap.find(iDofEdge);
						int iDof = icmp*numTotalDof2 + iter->second;
						// "Interpolating" data at the nodal point
						Real value = _vec1( icmp*numTotalDof1 + iDofEdge );
						// now what to what boundary node ( in the solid numerotation ) should we send this value ?
						//                                            std::cout << " transfering " << value << " from P2 " << edgeID << " to P2 " << iDof  << " ( " << iter->second << " ) " << std::endl;
						// Updating interpolated mesh velocity
						_vec2.checkAndSet( iDof, value );
					}
				}
			}
			else
			{
				for ( ID l = 1; l <= nDofPerEdge1; ++l )
				{
					// Loop on data vector components
					for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
					{
						//
						// ID of the 1st node of the edge
						ID node1 = _fespace1.mesh()->edgeList(iEdge).point(1).id();
						iter = locDofMap.find(node1);
						// ID of the 1st dof of the edge in the solid numerotation
						int iDof1 = icmp*numTotalDof2 + iter->second;
						value = 0.5*_vec1( icmp*numTotalDof1 + iDofEdge  ) + _vec2[iDof1];
						//                                            std::cout << " transfering " << value << " from P2 " << iDofEdge << " to P1 " << iDof1 << " ( " << iter->second << " ) " << std::endl;
						_vec2.checkAndSet( iDof1, value );
						//
						// ID of the 2nd node of the edge
						ID node2 = _fespace1.mesh()->edgeList(iEdge).point(2).id();
						iter = locDofMap.find(node2);
						// ID of the 2nd dof of the edge in the solid numerotation
						int iDof2 = icmp*numTotalDof2 + iter->second;
						value = 0.5*_vec1( icmp*numTotalDof1 + iDofEdge ) + _vec2[iDof2];
						//                                            std::cout << " transfering " << value << " from P2 " << iDofEdge << " to P1 " << iDof2 << " ( " << iter->second << " ) " << std::endl;
						_vec2.checkAndSet( iDof2, value );
						// now what to what boundary node ( in the solid numerotation ) should we send this value ?
						//std::cout << "" << std::endl;
						// Updating interpolated mesh velocity
						//                                                     _vec2.checkAndSet( iDof2, value );
						//                                                     std::cout << std::endl;
					}
				}
			}
		}
	}
    else if (nDofPerEdge2)
	{
		// The first FESpace has no dofs on the edge
		// The second FESpace has dofs on the edge
		// We need to interpolate the vertex dofs on the edge dofs

		// First we need to go through the FS interface boundary edges on the second mesh

		// loop on boundary edges
		for ( ID iEdge = 1; iEdge <= nBEdges2; ++iEdge )
		{
			if ((int)_fespace2.mesh()->edgeList(iEdge).marker() != M_fluidInterfaceFlag) continue;
			// Now that we have an edge on the FS interface, let's get its ID
			ID edgeID = _fespace2.mesh()->edgeList(iEdge).id();
			// dof position of the edge since unknowns are store using [ node | edges | faces | volumes ]
			int iDofEdge = numTotalVert2 + edgeID;

			for ( ID l = 1; l <= nDofPerEdge1; ++l )
			{
				// Loop on data vector components
				for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
				{
					// interpolation of the nodal values in the middle of the segement
					ID node1 = _fespace2.mesh()->edgeList(iEdge).point(1).id();
					ID node2 = _fespace2.mesh()->edgeList(iEdge).point(2).id();
					value = 0.5*(_vec2(icmp*numTotalDof2 + node1) + _vec2(icmp*numTotalDof2 + node2));
					_vec2.checkAndSet(iDofEdge, value);
				}
			}
		}
	}


    //         // Face based Dof
//         if ( nDofpF )
//         {

//             // loop on element faces
//             for ( ID iFa = 1; iFa <= nElemF; ++iFa )
//             {

//                 // Loop on number of Dof per face
//                 for ( ID l = 1; l <= nDofpF; ++l )
//                 {

//                     lDof = nDofElemE + nDofElemV + ( iFa - 1 ) * nDofpF + l; // Local dof in the adjacent Element

//                     // Nodal coordinates
//                     x = _fespace1.refFE().xi( lDof - 1 );
//                     y = _fespace1.refFE().eta( lDof - 1 );
//                     z = _fespace1.refFE().zeta( lDof - 1 );

//                     // Loop on data vector components
//                     for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
//                     {

//                         // Interpolating data at the nodal point
//                         Real __sum = 0;
//                         for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local Dof on the adjacent element
//                             __sum += wLoc( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi( idof, x, y, z ); // Problem here with P2

//                         // Updating interpolating vector
//                         int iDof = icmp * _fespace1.dof().numTotalDof() + _fespace1.dof().localToGlobal( iElem, lDof + 1);
//                         _vec2.checkAndSet( iDof ,__sum);
//                     }
//                 }
//             }
//         }

//         // Element based Dof
//         // Loop on number of Dof per Element
//         for ( ID l = 1; l <= nDofpEl; ++l )
//         {
//             lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

//             // Nodal coordinates
//             x = _fespace1.refFE().xi( lDof - 1 );
//             y = _fespace1.refFE().eta( lDof - 1 );
//             z = _fespace1.refFE().zeta( lDof - 1 );

//             // Loop on data vector components
//             for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
//             {

//                 // Interpolating data at the nodal point
//                 Real __sum = 0;
//                 for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local Dof on the adjacent element
//                     __sum += wLoc( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi( idof, x, y, z );

//                 // Updating interpolating vector

//                 int iDof = icmp * _fespace1.dof().numTotalDof() + _fespace1.dof().localToGlobal( elemId, lDof );
//                 _vec2.checkAndSet( iDof, __sum);
//             }
//         }
}

} // Namespace LifeV

#endif /* TWODIM */
