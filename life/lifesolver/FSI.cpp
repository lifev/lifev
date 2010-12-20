//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER
//  \include ../../../mathcard/testsuite/test_monolithic/fluidstructure.dox


/*!
 \include ../../doc/api/bibliography/newton
 \include ../../doc/api/bibliography/fluidstructure

    @file
    @brief Pure virtual operator class for FSI solvers

    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gilles Fourestey <fourestey@cscs.ch>
    @contributor Paolo Crosetto <paolo.crosetto@epfl.ch>
    @maintainer  Paolo Crosetto <paolo.crosetto@epfl.ch>

    @date 10-12-2010

    This is the base class for the FSI solvers in LifeV. It contains the methods to evaluate the residual and compute the
    Jacobian matrix, which make it suited for the generalized Newton method implemente in NonlinearRichardson. The fluid
    and structure classes are member of this class and different formulations (e.g. Monolithic \ref CDFQ , segregated
    Newton \ref FM05 , Dirichlet--Neumann \ref DDFQ06 , Robin Neumann \ref BNV08 )

 */

#include <life/lifecore/life.hpp>
#include <life/lifesolver/FSI.hpp>

namespace LifeV
{

// ===================================================
//! Constructors & Destructors
// ===================================================
FSI::FSI():
    M_mesh                               ( ),
    M_uFESpace                           ( ),
    M_pFESpace                           ( ),
    M_dFESpace                           ( ),
    M_mmFESpace                          ( ),
    M_fluidMesh                          ( ),
    M_solidMesh                          ( ),
    M_fluidMeshPart                      ( ),
    M_solidMeshPart                      ( ),
    M_BCh_u                              ( ),
    M_BCh_d                              ( ),
    M_BCh_mesh                           ( ),
    M_BCh_du                             ( ),
    M_BCh_du_inv                         ( ),
    M_BCh_dz                             ( ),
    M_BCh_dz_inv                         ( ),
    M_BCh_dp                             ( ),
    M_BCh_dp_inv                         ( ),
    M_fluid                              ( ),
    M_solid                              ( ),
    M_meshMotion                         ( ),
    //     M_fluidLin                           ( ),
    //     M_solidLin                           ( ),
    M_bdf                                ( ),
    M_dataFile                           ( ),
    M_dataMeshFluid                      ( new DataMesh()),
    M_dataMeshSolid                      ( new DataMesh()),
    M_data                               ( ),
    M_fluidInterfaceMap                  ( ),
    M_solidInterfaceMap                  ( ),
    M_fluidInterfaceMapOnZero            ( ),
    M_solidInterfaceMapOnZero            ( ),
    M_dofFluidToStructure                ( ),
    M_dofStructureToFluid                ( ),
    M_dofStructureToSolid                ( ),
    M_dofStructureToHarmonicExtension    ( ),
    M_dofHarmonicExtensionToFluid        ( ),
    M_dofFluid                           ( ),
    M_dofSolid                           ( ),
    M_dofFluidInv                        ( ),
    M_dofSolidInv                        ( ),
    M_bcvFluidInterfaceDisp              ( ),
    M_bcvFluidLoadToStructure            ( ),
    M_bcvSolidLoadToStructure            ( ),
    M_bcvStructureToFluid                ( ),
    M_bcvStructureDispToFluid            ( ),
    M_bcvStructureDispToSolid            ( ),
    M_bcvStructureDispToHarmonicExtension( ),
    M_bcvHarmonicExtensionVelToFluid     ( ),
    M_bcvDerHarmonicExtensionVelToFluid  ( ),
    M_bcvDerFluidLoadToStructure         ( ),
    M_bcvDerFluidLoadToFluid             ( ),
    M_bcvDerStructureDispToSolid         ( ),
    M_bcfRobinOuterWall                  ( ),
    M_lambdaFluid                        ( ),
    M_lambdaFluidRepeated                ( ),
    M_lambda                             ( ),
    M_lambdaDot                          ( ),
    M_un                                 ( ),
    M_rhs                                ( ),
    M_Alphaf                             ( ), //vector_Type, for alphaf robin
    M_AlphafCoef                         ( 0 ),
    M_betamedio                          ( ),
    M_fluxes                             ( 0 ),
    M_epetraComm                         ( ),
    M_epetraWorldComm                    ( ),
    //begin of private members
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
    M_mpi                                ( true ),
    M_isFluid                            ( false ),
    M_isSolid                            ( false ),
    M_linearFluid                        ( false ),
    M_linearSolid                        ( false ),
    M_fluidLeader                        ( ),
    M_solidLeader                        ( )
{
}

FSI::~FSI()
{
}


// ===================================================
//! Virtual Public Methods
// ===================================================

void
FSI::setDataFile( const GetPot& dataFile )
{

    M_fluidMesh.reset(new mesh_Type());
    M_dataMeshFluid->setup(dataFile, "fluid/space_discretization");
    readMesh(*M_fluidMesh, *M_dataMeshFluid);

    M_solidMesh.reset(new mesh_Type());
    M_dataMeshSolid->setup(dataFile, "solid/space_discretization");
    readMesh(*M_solidMesh, *M_dataMeshSolid);

    M_dataFile = dataFile;
}

void
FSI::setupFEspace()
{
    Displayer disp(M_epetraComm);
    disp.leaderPrint("FSI-  Setting RefFE and QuadRule ...           ");

    std::string uOrder = M_data->dataFluid()->uOrder();
    std::string pOrder = M_data->dataFluid()->pOrder();
    std::string dOrder = M_data->dataSolid()->getOrder();

    const RefFE*    refFE_vel(0);
    const QuadRule* qR_vel(0);
    const QuadRule* bdQr_vel(0);

    const RefFE*    refFE_press(0);
    const QuadRule* qR_press(0);
    const QuadRule* bdQr_press(0);

    const RefFE*    refFE_struct(0);
    const QuadRule* qR_struct(0);
    const QuadRule* bdQr_struct(0);

    if ( uOrder.compare("P2") == 0 )
    {
        refFE_vel = &feTetraP2;
        qR_vel    = &quadRuleTetra15pt; // DoE 5
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( uOrder.compare("P1") == 0 )
    {
        refFE_vel = &feTetraP1;
        qR_vel    = &quadRuleTetra4pt;  // DoE 2
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( uOrder.compare("P1Bubble") == 0 )
    {
        refFE_vel = &feTetraP1bubble;
        qR_vel    = &quadRuleTetra64pt;  // DoE 2
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else
    {
        ERROR_MSG(uOrder + " velocity FE not implemented yet.");
    }

    if ( pOrder.compare("P2") == 0 )
    {
        refFE_press = &feTetraP2;
        qR_press    = qR_vel; // DoE 5
        bdQr_press  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( pOrder.compare("P1") == 0 )
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

    if ( dOrder.compare("P2") == 0 )
    {
        refFE_struct = &feTetraP2;
        qR_struct    = &quadRuleTetra15pt; // DoE 5
        bdQr_struct  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( dOrder.compare("P1") == 0 )
    {
        refFE_struct = &feTetraP1;
        qR_struct    = &quadRuleTetra4pt;  // DoE 2
        bdQr_struct  = &quadRuleTria3pt;   // DoE 2
    }
    else
    {
        ERROR_MSG(dOrder + " structure FE not implemented yet.");
    }
    disp.leaderPrint("done\n");

    disp.leaderPrint("FSI-  Building fluid FESpace ...               \n");
    if (this->isFluid())
    {

        M_mmFESpace.reset(new FESpace<mesh_Type, EpetraMap>(*M_fluidMeshPart,
                                                            //dOrder,
                                                            *refFE_struct,
                                                            *qR_struct,
                                                            *bdQr_struct,
                                                            3,
                                                            M_epetraComm));

        M_uFESpace.reset( new FESpace<mesh_Type, EpetraMap>(*M_fluidMeshPart,
                                                            //uOrder,
                                                            *refFE_vel,
                                                            *qR_vel,
                                                            *bdQr_vel,
                                                            3,
                                                            M_epetraComm));

        M_pFESpace.reset( new FESpace<mesh_Type, EpetraMap>(*M_fluidMeshPart,
                                                            //pOrder,
                                                            *refFE_press,
                                                            *qR_press,
                                                            *bdQr_press,
                                                            1,
                                                            M_epetraComm));
    }
    else
    {
        M_mmFESpace.reset(new FESpace<mesh_Type, EpetraMap>(M_fluidMesh,
                                                            //dOrder,
                                                            *refFE_struct,
                                                            *qR_struct,
                                                            *bdQr_struct,
                                                            3,
                                                            M_epetraComm));

        M_uFESpace.reset( new FESpace<mesh_Type, EpetraMap>(M_fluidMesh,
                                                            //uOrder,
                                                            *refFE_vel,
                                                            *qR_vel,
                                                            *bdQr_vel,
                                                            3,
                                                            M_epetraComm));

        M_pFESpace.reset( new FESpace<mesh_Type, EpetraMap>(M_fluidMesh,
                                                            //pOrder,
                                                            *refFE_press,
                                                            *qR_press,
                                                            *bdQr_press,
                                                            1,
                                                            M_epetraComm));
    }
    M_epetraWorldComm->Barrier();

    disp.leaderPrint("FSI-  Building solid FESpace ...               \n");
    if (this->isSolid())
    {
        M_dFESpace.reset( new FESpace<mesh_Type, EpetraMap>( *M_solidMeshPart,
                                                             dOrder,
                                                             //*refFE_struct,
                                                             //*qR_struct,
                                                             //*bdQr_struct,
                                                             3,
                                                             M_epetraComm));
    }
    else
    {
        M_dFESpace.reset(new FESpace<mesh_Type, EpetraMap>(M_solidMesh,
                                                           //dOrder,
                                                           *refFE_struct,
                                                           *qR_struct,
                                                           *bdQr_struct,
                                                           3,
                                                           M_epetraComm));
    }
    M_epetraWorldComm->Barrier();
}


void
FSI::partitionMeshes()
{
    if (this->isFluid())
    {
        M_fluidMeshPart.reset(new  MeshPartitioner< mesh_Type > (M_fluidMesh, M_epetraComm));
    }
    if (this->isSolid())
    {
        M_solidMeshPart.reset( new  MeshPartitioner< mesh_Type > ( M_solidMesh, M_epetraComm ) );
    }

}

#ifdef HAVE_HDF5
void
FSI::partitionMeshes( meshFilter_Type& fluidMeshFilter, meshFilter_Type& solidMeshFilter )
{
    M_fluidMesh = fluidMeshFilter.getMeshPartition();
    M_solidMesh = solidMeshFilter.getMeshPartition();
}
#endif


void
FSI::setupDOF( void )
{
    Displayer disp(M_epetraWorldComm);
    disp.leaderPrint("FSI: setting DOF ... " );
    Dof uDof(*M_fluidMesh, M_uFESpace->refFE()); // velocity dof related to unpartitioned mesh
    Dof dDof(*M_solidMesh, M_dFESpace->refFE()); // velocity dof related to unpartitioned mesh

    M_dofFluidToStructure                .reset( new DofInterface3Dto3D );
    M_dofStructureToFluid                .reset( new DofInterface3Dto3D );
    M_dofStructureToSolid                .reset( new DofInterface3Dto3D );
    M_dofStructureToHarmonicExtension    .reset( new DofInterface3Dto3D );
    M_dofHarmonicExtensionToFluid        .reset( new DofInterface3Dto3D );
    M_dofFluid                           .reset( new DofInterface3Dto2D );
    M_dofSolid                           .reset( new DofInterface3Dto2D );
    M_dofFluidInv                        .reset( new DofInterface3Dto2D );
    M_dofSolidInv                        .reset( new DofInterface3Dto2D );
    M_bcvFluidInterfaceDisp              .reset( new  BCVectorInterface );
    M_bcvFluidLoadToStructure            .reset( new  BCVectorInterface );
    M_bcvSolidLoadToStructure            .reset( new  BCVectorInterface );
    M_bcvStructureToFluid                .reset( new  BCVectorInterface );
    M_bcvStructureDispToFluid            .reset( new  BCVectorInterface );
    M_bcvStructureDispToSolid            .reset( new  BCVectorInterface );
    M_bcvStructureDispToHarmonicExtension.reset( new  BCVectorInterface );
    M_bcvHarmonicExtensionVelToFluid     .reset( new  BCVectorInterface );
    M_bcvDerHarmonicExtensionVelToFluid  .reset( new  BCVectorInterface );
    M_bcvDerFluidLoadToStructure         .reset( new  BCVectorInterface );
    M_bcvDerFluidLoadToFluid             .reset( new  BCVectorInterface );
    M_bcvDerStructureDispToSolid         .reset( new  BCVectorInterface );

    M_dofFluidToStructure->setup(   M_dFESpace->refFE(), M_dFESpace->dof(),
                                    M_uFESpace->refFE(), uDof );
    M_dofFluidToStructure->update( *M_dFESpace->mesh(),  M_data->structureInterfaceFlag(),
                                   *M_fluidMesh,         M_data->fluidInterfaceFlag(),
                                   M_data->interfaceTolerance(),
                                   M_data->structureInterfaceVertexFlag() );

    //here the solid mesh must be non partitioned in the monolithic case
    M_dofStructureToHarmonicExtension->setup(   M_uFESpace->refFE(), M_uFESpace->dof(),
                                                M_dFESpace->refFE(), dDof );
    M_dofStructureToHarmonicExtension->update( *M_uFESpace->mesh(),  M_data->fluidInterfaceFlag(),
                                               *M_solidMesh, M_data->structureInterfaceFlag(),
                                               M_data->interfaceTolerance(),
                                               M_data->fluidInterfaceVertexFlag() );

    M_dofStructureToSolid->setup(   M_dFESpace->refFE(), M_dFESpace->dof(),
                                    M_dFESpace->refFE(), M_dFESpace->dof() );
    M_dofStructureToSolid->update( *M_dFESpace->mesh(),  M_data->structureInterfaceFlag(),
                                   *M_dFESpace->mesh(),  M_data->structureInterfaceFlag(),
                                   M_data->interfaceTolerance(),
                                   M_data->structureInterfaceVertexFlag() );

    M_dofStructureToFluid->setup(   M_uFESpace->refFE(), M_uFESpace->dof(), //modifica matteo FSI
                                    M_dFESpace->refFE(), dDof );
    M_dofStructureToFluid->update( *M_uFESpace->mesh(), M_data->fluidInterfaceFlag(),
                                   *M_solidMesh,        M_data->structureInterfaceFlag(),
                                   M_data->interfaceTolerance(),
                                   M_data->fluidInterfaceVertexFlag());

    M_dofHarmonicExtensionToFluid->setup(   M_uFESpace->refFE(),  M_uFESpace->dof(),
                                            M_uFESpace->refFE(),  M_uFESpace->dof() );
    M_dofHarmonicExtensionToFluid->update( *M_uFESpace->mesh(),  M_data->fluidInterfaceFlag(),
                                           *M_uFESpace->mesh(),  M_data->fluidInterfaceFlag(),
                                           M_data->interfaceTolerance(),
                                           M_data->fluidInterfaceVertexFlag());

    M_epetraWorldComm->Barrier();
    disp.leaderPrint("done\n");

    createInterfaceMaps( M_dofStructureToHarmonicExtension->localDofMap() );
}

void
FSI::setupFluidSolid( void )
{
    setupFluidSolid(imposeFlux());
}

void
FSI::setupFluidSolid( UInt const fluxes )
{
    if ( this->isFluid() )
    {
        //M_fluxes = imposeFlux();

        M_meshMotion.reset( new meshMotion_Type(               *M_mmFESpace,             M_epetraComm ) );
        M_fluid.reset(      new fluid_Type(      M_data->dataFluid(), *M_uFESpace, *M_pFESpace, M_epetraComm, fluxes ) );
        M_solid.reset( solid_Type::StructureSolverFactory::instance().createObject( M_data->dataSolid()->getSolidType( ) ) );
        M_solid->setup( M_data->dataSolid(), M_dFESpace, M_epetraComm );

//         if ( M_linearFluid )
//             M_fluidLin.reset( new FSI::fluidlin_raw_type( *M_data->dataFluid(), *M_uFESpace, *M_pFESpace, M_epetraComm ) );

//         if ( M_linearSolid )
//             M_solidLin.reset( new FSI::solidlin_raw_type( *M_data->dataSolid(), *M_dFESpace, M_epetraComm ) );

        //Vector initialization
        M_rhs.reset( new vector_Type( M_fluid->getMap() ) );
    }

    if ( this->isSolid() )
    {
//         M_fluid.reset( new fluid_Type( *M_data->dataFluid(), *M_uFESpace, *M_pFESpace, M_epetraComm ) );
        M_meshMotion.reset( new meshMotion_Type(               *M_mmFESpace, M_epetraComm ) );
        M_solid.reset(solid_Type::StructureSolverFactory::instance().createObject( M_data->dataSolid()->getSolidType( ) ) );
        M_solid->setup( M_data->dataSolid(), M_dFESpace,  M_epetraComm );

//         if ( M_linearFluid )
//             M_fluidLin.reset( new FSI::fluidlin_raw_type( *M_data->dataFluid(), *M_uFESpace, *M_pFESpace, M_epetraComm ) );

//         if ( M_linearSolid )
//             M_solidLin.reset( new FSI::solidlin_raw_type( *M_data->dataSolid(), *M_dFESpace, M_epetraComm ) );
    }

    M_epetraWorldComm->Barrier();
}

void
FSI::setupSystem( void )
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
        M_solid->setDataFromGetPot( M_dataFile );
//         if (M_linearSolid)
//             M_solidLin->setUp( dataFile );
    }

    M_epetraWorldComm->Barrier();
}


void
FSI::buildSystem()
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
FSI::updateSystem( )
{
    shiftSolution();

    if ( this->isFluid() )
    {
        M_meshMotion->updateSystem();

        transferMeshMotionOnFluid(M_meshMotion->disp(), *this->M_dispFluidMeshOld);

        if ( M_fluid->solution().get() )
            M_un.reset( new vector_Type( *M_fluid->solution() ) );
	M_bdf->updateRHSContribution( M_data->dataFluid()->dataTime()->timeStep() );
        *M_rhs = M_fluid->matrixMass()*M_bdf->rhsContributionFirstDerivative();
    }

    if ( this->isSolid() )
    {
        this->M_solid->updateSystem();
    }
   this->couplingVariableExtrap( );
}



void FSI::couplingVariableExtrap( )
{
    *M_lambda      = lambdaSolid();
    if (!M_lambdaDot.get())
    {
        M_lambdaDot.reset        ( new vector_Type(*M_fluidInterfaceMap, Unique) );
        *M_lambda     += M_data->dataFluid()->dataTime()->timeStep()*lambdaDotSolid();
    }
    else
    {
        *M_lambda     += 1.5*M_data->dataFluid()->dataTime()->timeStep()*lambdaDotSolid(); // *1.5
        *M_lambda     -= M_data->dataFluid()->dataTime()->timeStep()*0.5*(*M_lambdaDot);
    }

    *M_lambdaDot   = lambdaDotSolid();
    displayer().leaderPrint("FSI-  norm( disp  ) init =                     ", M_lambda->normInf(), "\n" );
    displayer().leaderPrint("FSI-  norm( velo )  init =                     ", M_lambdaDot->normInf(), "\n");
}

void
FSI::shiftSolution()
{
    if ( this->isFluid() )
    {
        this->M_bdf->shiftRight( *M_fluid->solution() );
    }
}

UInt
FSI::imposeFlux( void )
{
    if ( this->isFluid() )
    {
        std::vector<bcName_Type> fluxVector = M_BCh_u->findAllBCWithType( Flux );
        UInt numLM = static_cast<UInt>( fluxVector.size() );

        UInt offset = M_uFESpace->map().map(Unique)->NumGlobalElements()
                      + M_pFESpace->map().map(Unique)->NumGlobalElements();

        for ( UInt i = 0; i < numLM; ++i )
            M_BCh_u->setOffset( fluxVector[i], offset + i );
        return numLM;
    }
    else
        return 0;
}

// ===================================================
//! Public Methods
// ===================================================


void FSI::initializeBDF( const vector_Type& un )
{
  M_bdf.reset( new BdfT<vector_Type>( ));
  M_bdf->setup(M_data->dataFluid()->dataTime()->orderBDF() ) ;
 M_bdf->setInitialCondition( un );
}

void FSI::createInterfaceMaps( std::map<ID, ID> const& locDofMap )
{
    Displayer disp(M_epetraWorldComm);
    disp.leaderPrint("FSI-  Building fluid variables ...             ");

    // now we build the sigma and lambda variables on each proc


    std::vector<int> dofInterfaceFluid;

    //is the interface map between HE (first) and solid (second)
    if ( this->isFluid() )
    {
        //std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
        dofInterfaceFluid.reserve( locDofMap.size() );

        for (UInt dim = 0; dim < nDimensions; ++dim)
            for ( iterator_Type i = locDofMap.begin(); i != locDofMap.end(); ++i )
                dofInterfaceFluid.push_back(i->second + dim * M_dFESpace->dof().numTotalDof()); // in solid numerotation
    }

    int* pointerToDofs(0);
    if (dofInterfaceFluid.size() > 0)
        pointerToDofs = &dofInterfaceFluid[0];

    M_fluidInterfaceMap.reset( new EpetraMap( -1,
                                              static_cast<int>(dofInterfaceFluid.size()),
                                              pointerToDofs,
                                              1,
                                              M_epetraWorldComm ));
    disp.leaderPrint("done\n");
    M_epetraWorldComm->Barrier();

    disp.leaderPrint("FSI-  Building solid variables ...             ");

    std::vector<int> dofInterfaceSolid;

    if (this->isSolid())
    {
        //std::map<ID, ID> const& locDofMap = M_dofFluidToStructure->locDofMap();
        dofInterfaceSolid.reserve( locDofMap.size() );
        //std::cout << "solid" << std::endl;
        for (UInt dim = 0; dim < nDimensions; ++dim)
            for ( iterator_Type i = locDofMap.begin(); i != locDofMap.end(); ++i )
            {
                dofInterfaceSolid.push_back(i->second/*Simone: first*/ + dim * M_dFESpace->dof().numTotalDof()); // in solid numerotation
            }
    }


    pointerToDofs = 0;
    if (dofInterfaceSolid.size() > 0)
        pointerToDofs = &dofInterfaceSolid[0];

    M_solidInterfaceMap.reset( new EpetraMap( -1,
                                              static_cast<int>(dofInterfaceSolid.size()),
                                              pointerToDofs,
                                              1,
                                              M_epetraWorldComm ));

    M_epetraWorldComm->Barrier();
    disp.leaderPrint("done\n");

    disp.leaderPrint("FSI-  Variables initialization ...             \n");

    //variablesInit( refFE_struct, bdQr_struct, qR_struct);
    variablesInit( M_data->dataSolid()->getOrder() );

    M_epetraWorldComm->Barrier();
}



void
FSI::initializeFluid( const vector_Type& velAndPressure,
                              const vector_Type& displacement )
{
    this->fluid().initialize( velAndPressure );
    this->meshMotion().initialize( displacement );
    this->moveMesh( displacement);
}



void
FSI::initializeSolid( vectorPtr_Type displacement,
                              vectorPtr_Type velocity )
{
    this->solid().initialize( displacement, velocity);
}



void
FSI::moveMesh( const vector_Type& dep )
{
    displayer().leaderPrint("FSI-  Moving the mesh ...                      ");
    M_fluidMeshPart->meshPartition()->moveMesh(dep,  this->M_mmFESpace->dof().numTotalDof());
    displayer().leaderPrint( "done\n" );
    M_fluid->setRecomputeMatrix( true );
}



void
FSI::transferFluidOnInterface(const vector_Type &_vec1, vector_Type &_vec2)
{
    // e.g.: vec1=M_fluid->residual(), vec2=M_sigmaFluid
//     _vec2 = ZeroVector(_vec2.size());

    if (_vec1.mapType() == Unique)
    {
        vector_Type const  vec1Repeated(_vec1, Repeated);
        transferFluidOnInterface(vec1Repeated, _vec2);
        return;
    }

    if (_vec2.mapType() == Repeated)
    {
        vector_Type  vec2Unique(_vec2, Unique);
        transferFluidOnInterface(_vec1, vec2Unique);
        _vec2 = vec2Unique;
        return;
    }

    _vec2 *= 0;

    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->localDofMap();

    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();
    int numTotalDofFluid = M_uFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::const_iterator iterator_Type;

    for (UInt dim = 0; dim < nDimensions; ++dim)
        for ( iterator_Type it = locDofMap.begin(); it != locDofMap.end(); ++it )
        {
//                 std::cout <<  it->second + dim*numTotalDofFluid << " to "
//                           <<  it->first + dim*numTotalDofSolid << " : "
//                           << _vec1[it->second + dim*numTotalDofFluid] << std::endl;
            _vec2.setCoefficient( it->second + dim*numTotalDofSolid,
                                  _vec1[it->first + dim*numTotalDofFluid] );
        }
}



//works in serial but no yet in parallel
void
FSI::transferSolidOnFluid(const vector_Type &_vec1, vector_Type &_vec2)
{
    //    e.g.: vec2=M_fluid->residual(), vec1=M_sigmaFluid
    //     _vec2 = ZeroVector(_vec2.size());

    if (_vec1.mapType() == Unique)
    {
        vector_Type const  vec1Repeated(_vec1, Repeated);
        transferSolidOnInterface(vec1Repeated, _vec2);
        return;
    }

    if (_vec2.mapType() == Repeated)
    {
        vector_Type  vec2Unique(_vec2, Unique);
        transferSolidOnInterface(_vec1, vec2Unique);
        _vec2 = vec2Unique;
        return;
    }

    _vec2 *= 0;

    std::map<ID, ID> const& locDofMap = M_dofFluidToStructure->localDofMap();


    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();
    int numTotalDofFluid = M_uFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::const_iterator iterator_Type;

    for (UInt dim = 0; dim < nDimensions; ++dim)
        for ( iterator_Type it = locDofMap.begin(); it != locDofMap.end(); ++it )
        {
            _vec2.setCoefficient( it->second + dim*numTotalDofFluid,
                                  _vec1[it->first + dim*numTotalDofSolid]);
        }
}



void
FSI::transferSolidOnInterface(const vector_Type &_vec1, vector_Type &_vec2)
{
    /* e.g.:
       vec1 (Unique)       vec2 (Repeated)  (On different EpetraMaps) (Changing now to Unique on both)
       M_solid->disp()     M_lambdaSolid
       M_solid->vel()      M_lambdaDotSolid
       M_solid->residual() M_sigmaSolid
    */

    if (_vec1.mapType() == Unique)
    {
        vector_Type const  vec1Repeated(_vec1, Repeated);
        transferSolidOnInterface(vec1Repeated, _vec2);
        return;
    }

    if (_vec2.mapType() == Repeated)
    {
        vector_Type  vec2Unique(_vec2, Unique);
        transferSolidOnInterface(_vec1, vec2Unique);
        _vec2 = vec2Unique;
        return;
    }

    _vec2 *= 0;

    std::map<ID, ID> const& locDofMap = M_dofStructureToSolid->localDofMap();

    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::const_iterator iterator_Type;

    for (UInt dim = 0; dim < nDimensions; ++dim)
        for ( iterator_Type it = locDofMap.begin(); it != locDofMap.end(); ++it )
        {
            _vec2.setCoefficient( it->second + dim*numTotalDofSolid,
                                  _vec1[it->first + dim*numTotalDofSolid] );
        }
}

void
FSI::transferInterfaceOnSolid(const vector_Type& _vec1, vector_Type& _vec2)
{
    /* e.g.:
       vec2                vec1
       M_solid->disp()     M_lambdaSolid
       M_solid->vel()      M_lambdaDotSolid
       M_solid->residual() M_sigmaSolid
    */

    if (_vec1.mapType() == Unique)
    {
        vector_Type const  vec1Repeated(_vec1, Repeated);
        transferInterfaceOnSolid(vec1Repeated, _vec2);
        return;
    }

    if (_vec2.mapType() == Repeated)
    {
        vector_Type  vec2Unique(_vec2, Unique);
        transferInterfaceOnSolid(_vec1, vec2Unique);
        _vec2 = vec2Unique;
        return;
    }

    _vec2 *= 0;

    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->localDofMap();
    int numTotalDofSolid = M_dFESpace->dof().numTotalDof();

    typedef std::map<ID, ID>::const_iterator iterator_Type;

    for (UInt dim = 0; dim < nDimensions; ++dim)
        for ( iterator_Type it = locDofMap.begin(); it != locDofMap.end(); ++it )
        {
            _vec2.setCoefficient( it->second + dim*numTotalDofSolid,
                                  _vec1[it->second + dim*numTotalDofSolid] );
        }
}

void
FSI::bcManageVectorRHS( const fluidBchandlerPtr_Type& bch, vector_Type& rhs )
{
    if ( !bch->bcUpdateDone() || M_fluid->recomputeMatrix() )
        bch->bcUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );

    bcManageVector( rhs, *M_uFESpace->mesh(), M_uFESpace->dof(),  *bch, M_uFESpace->feBd(), 1., 0. );
}

void
FSI::setAlphafCoef( )
{
    Real h=0.1, R=0.5;
    M_AlphafCoef  = 2*(this->dataSolid()->getRho()*h)/this->dataFluid()->dataTime()->timeStep();
    M_AlphafCoef += h*this->dataSolid()->getYoung(0)*this->dataFluid()->dataTime()->timeStep() /
                    (2*pow(R,2) *(1-pow(dataSolid()->getPoisson(0),2)));
}

void
FSI::setStructureToFluidParametres()
{
    this->setAlphafCoef();
    this->setAlphaf();

    if (M_Alphaf.get()==0)
    {
        this->setAlphafCoef();
        M_bcvStructureToFluid->setRobinCoeff(M_AlphafCoef);
        M_bcvStructureToFluid->setBetaCoeff(M_AlphafCoef);
    }
    else
    {
        M_bcvStructureToFluid->setRobinCoeffVector(this->Alphaf());
        M_bcvStructureToFluid->setBetaCoeffVector(this->Alphaf());
    }
}

// ===================================================
//! Display Methods
// ===================================================
bool
FSI::isLeader() const
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
FSI::displayer()
{
    if ( isFluid() &&  M_fluid.get())
        return M_fluid->getDisplayer();

    if ( !isSolid() || !M_solid.get() )
        std::cout << "ERROR: displayer not ready" << std::endl;

    return M_solid->getDisplayer();
}





// ===================================================
//! Get Functions
// ===================================================
/*
FSI::vector_Type
FSI::displacementOnInterface()
{

//     vector_Type dispOnInterface();
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
FSI::setComm( const commPtr_Type& comm,
                      const commPtr_Type& worldComm )
{
    M_epetraComm       = comm;
    M_epetraWorldComm  = worldComm;
}


void
FSI::setFluid( const fluidPtr_Type& fluid, const meshMotionPtr_Type& meshmotion )
{
    M_fluid = fluid;
    M_meshMotion = meshmotion;
    M_isFluid = true;
}



void
FSI::setSolid( const solidPtr_Type& solid )
{
    M_solid = solid;
    M_isSolid = true;
}



void
FSI::setFluidBC( const fluidBchandlerPtr_Type& bc_fluid )
{
    if ( isFluid() )
    {
        M_BCh_u = bc_fluid;
    }
}



void
FSI::setHarmonicExtensionBC( const fluidBchandlerPtr_Type& bc_he )
{
    if ( isFluid() )
    {
        M_BCh_mesh = bc_he;
        //M_meshMotion->setBC(*M_BCh_mesh);
    }
}



void
FSI::setSolidBC( const solidBchandlerPtr_Type& bc_solid )
{
    if ( isSolid() )
    {
        M_BCh_d = bc_solid;
    }
}



void
FSI::setLambdaFluid( const vector_Type& lambda )
{
    if ( lambda.mapType() == Unique )
        *M_lambdaFluid = lambda;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_lambdaFluidRepeated = *M_lambdaFluid;
}



void
FSI::setLambdaSolid( const vector_Type& lambda )
{
    if ( lambda.mapType() == Unique )
        *M_lambdaSolid = lambda;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_lambdaSolidRepeated = *M_lambdaSolid;
}



void
FSI::setLambdaSolidOld( const vector_Type& lambda )
{
    if ( lambda.mapType() == Unique )
        *M_lambdaSolidOld = lambda;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry
}



void
FSI::setLambdaDotSolid( const vector_Type& lambda )
{
    if ( lambda.mapType() == Unique )
        *M_lambdaDotSolid = lambda;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_lambdaDotSolidRepeated = *M_lambdaDotSolid;
}



void
FSI::setSigmaSolid( const vector_Type& sigma )
{
    if ( sigma.mapType() == Unique )
        *M_sigmaSolid = sigma;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_sigmaSolidRepeated = *M_sigmaSolid;
}


void
FSI::setSigmaFluid( const vector_Type& sigma )
{
    if ( sigma.mapType() == Unique )
        *M_sigmaFluid = sigma;
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_sigmaFluidRepeated = *M_sigmaFluid;

}



void
FSI::setMinusSigmaFluid( const vector_Type& sigma )
{
    if ( sigma.mapType() == Unique )
    {
        *M_minusSigmaFluid  = sigma;
        *M_minusSigmaFluid *= -1;
    }
    else // to be coded, I am not sure we need this functionality.
        assert(false); // if you get here, reformulate your problem in order to get a unique map as entry

    *M_minusSigmaFluidRepeated = *M_minusSigmaFluid;
}




void
FSI::setAlphafbcf( const bcFunction_Type& alphafbcf )
{
    vector_Type vec( M_fluid->velocityFESpace().map());
    M_fluid->velocityFESpace().interpolate(alphafbcf, vec, 0.0);
    *M_Alphaf = vec ;
}



void
FSI::setStructureDispToHarmonicExtension( const vector_Type& disp, UInt type )
{
    M_bcvStructureDispToHarmonicExtension->setup( disp,
                                                  M_dFESpace->dof().numTotalDof(),
                                                  M_dofStructureToHarmonicExtension,
                                                  type );
}



void
FSI::setStructureToFluid( const vector_Type& velo,  UInt type )
{
    M_bcvStructureToFluid->setup( velo,
                                  M_uFESpace->dof().numTotalDof(),
                                  M_dofHarmonicExtensionToFluid,
                                  type );
}



void
FSI::setStructureDispToFluid( const vector_Type& disp,  UInt type )
{
    M_bcvStructureDispToFluid->setup( disp,
                                      M_uFESpace->dof().numTotalDof(),
                                      M_dofStructureToFluid,
                                      type );
}



void
FSI::setStructureDispToSolid( const vector_Type& disp, UInt type )
{
    M_bcvStructureDispToSolid->setup( disp,
                                      M_dFESpace->dof().numTotalDof(),
                                      M_dofStructureToSolid,
                                      type );
}



void
FSI::setDerStructureDispToSolid( const vector_Type& ddisp, UInt type )
{
    M_bcvDerStructureDispToSolid->setup( ddisp,
                                         M_dFESpace->dof().numTotalDof(),
                                         M_dofStructureToSolid,
                                         type );
}



void
FSI::setSolidLoadToStructure( const vector_Type& load, UInt type )
{
    M_bcvSolidLoadToStructure->setup( load,
                                      M_dFESpace->dof().numTotalDof(),
                                      M_dofStructureToFluid,
                                      type );
}



void
FSI::setHarmonicExtensionVelToFluid( const vector_Type& vel, UInt type )
{
    M_bcvHarmonicExtensionVelToFluid->setup( vel,
                                             M_uFESpace->dof().numTotalDof(),
                                             M_dofHarmonicExtensionToFluid,
                                             type );
}



void
FSI::setDerHarmonicExtensionVelToFluid( const vector_Type& dvel, UInt type )
{
    M_bcvDerHarmonicExtensionVelToFluid->setup( dvel,
                                                M_uFESpace->dof().numTotalDof(),
                                                M_dofHarmonicExtensionToFluid,
                                                type );
}




void
FSI::setFluidLoadToStructure( const vector_Type& load, UInt type )
{
    M_bcvFluidLoadToStructure->setup( load,
                                      M_dFESpace->dof().numTotalDof(),
                                      M_dofStructureToSolid,
                                      type );
}



void
FSI::setDerFluidLoadToStructure( const vector_Type& dload, UInt type )
{
    M_bcvDerFluidLoadToStructure->setup( dload,
                                         M_dFESpace->dof().numTotalDof(),
                                         M_dofStructureToSolid,
                                         type );
}



void
FSI::setDerFluidLoadToFluid( const vector_Type& dload, UInt type )
{
    M_bcvDerFluidLoadToFluid->setup( dload,
                                     M_uFESpace->dof().numTotalDof(),
                                     M_dofHarmonicExtensionToFluid,
                                     type );
}

void FSI::setRobinOuterWall(function_Type const& dload, function_Type const& E)
{
    M_bcfRobinOuterWall.setFunctions_Robin(dload,
                                           E);
}


// ===================================================
//! Protected Methods
// ===================================================



void
FSI::variablesInit( const std::string& /*dOrder*/ )
//FSI::variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct)
{
    M_lambdaFluid.reset        ( new vector_Type(*M_fluidInterfaceMap, Unique) );
    M_lambda.reset             ( new vector_Type(*M_solidInterfaceMap, Unique) );
    M_lambdaDot.reset             ( new vector_Type(*M_solidInterfaceMap, Unique) );
    M_lambdaFluidRepeated.reset( new vector_Type(*M_fluidInterfaceMap, Repeated) );

    if ( this->isFluid() )
    {
        M_dispFluidMeshOld.reset( new vector_Type(M_uFESpace->map(), Repeated) );
        M_veloFluidMesh.reset   ( new vector_Type(M_uFESpace->map(), Repeated) );
        M_Alphaf.reset          ( new vector_Type(M_uFESpace->map(), Repeated));

        if ( M_linearFluid )
            M_derVeloFluidMesh.reset( new vector_Type(this->M_uFESpace->map(), Repeated) );
    }

    M_sigmaFluid.reset              ( new vector_Type(*M_fluidInterfaceMap, Unique) );
    M_sigmaFluidRepeated.reset      ( new vector_Type(*M_fluidInterfaceMap, Repeated) );
    M_minusSigmaFluid.reset         ( new vector_Type(*M_fluidInterfaceMap, Unique) );
    M_minusSigmaFluidRepeated.reset ( new vector_Type(*M_fluidInterfaceMap, Repeated) );

    M_lambdaSolid.reset   ( new vector_Type(*M_solidInterfaceMap, Unique) );
    M_lambdaSolidOld.reset( new vector_Type(*M_solidInterfaceMap, Unique) );
    M_lambdaDotSolid.reset( new vector_Type(*M_solidInterfaceMap, Unique) );
    M_sigmaSolid.reset    ( new vector_Type(*M_solidInterfaceMap, Unique) );

    M_lambdaSolidRepeated.reset   ( new vector_Type(*M_solidInterfaceMap, Repeated) );
    M_lambdaDotSolidRepeated.reset( new vector_Type(*M_solidInterfaceMap, Repeated) );
    M_sigmaSolidRepeated.reset    ( new vector_Type(*M_solidInterfaceMap, Repeated) );
}


void
FSI::transferMeshMotionOnFluid( const vector_Type& _vec1, vector_Type& _vec2 )
{
    //transferMeshMotionOnFluid should handle the repetition of the interface nodes.
    if (_vec1.mapType() == Unique)
    {
        vector_Type const  vec1Repeated(_vec1, Repeated);
        transferMeshMotionOnFluid(vec1Repeated, _vec2);
        return;
    }

    if (_vec2.mapType() == Repeated)
    {
        vector_Type  vec2Unique(_vec2, Unique);
        transferMeshMotionOnFluid(_vec1, vec2Unique);
        _vec2 = vec2Unique;
        return;
    }

    _vec2 *=0;

    interpolateVelocity(_vec1, _vec2);
    return;

}

void
FSI::interpolateVelocity( const vector_Type& _vec1, vector_Type& _vec2 )
{
    assert(_vec1.mapType() == Repeated);
    assert(_vec2.mapType() == Unique);

    typedef mesh_Type::VolumeShape GeoShape; // Element shape

    UInt nDofpV = M_uFESpace->refFE().nbDofPerVertex(); // number of Dof per vertex
    UInt nDofpE = M_uFESpace->refFE().nbDofPerEdge();   // number of Dof per edge
    UInt nDofpF = M_uFESpace->refFE().nbDofPerFace();   // number of Dof per face
    UInt nDofpEl = M_uFESpace->refFE().nbDofPerVolume(); // number of Dof per Volume

    UInt nElemV = GeoShape::S_numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::S_numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::S_numFaces;    // Number of element's faces

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
                        _vec2.setCoefficient( iDof ,__sum);

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
                        _vec2.setCoefficient( iDof ,__sum);

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
                        _vec2.setCoefficient( iDof ,__sum);
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
                _vec2.setCoefficient( iDof, __sum);
            }
        }
    }

}




// this will interpolate dofs values from fespace1 to fespace2
void
FSI::interpolateInterfaceDofs( const FESpace<mesh_Type, EpetraMap>& _fespace1,
                                       const vector_Type&                   _vec1,
                                       const FESpace<mesh_Type, EpetraMap>& _fespace2,
                                       vector_Type&                         _vec2,
                                       dofInterface3DPtr_Type&                _dofInterface)
{
    assert(_vec1.mapType() == Repeated);
    assert(_vec2.mapType() == Unique);

    typedef mesh_Type::VolumeShape GeoShape; // Element shape


    UInt nDofPerVert1  = _fespace1.refFE().nbDofPerVertex(); // number of Dof per vertex
    UInt nDofPerEdge1  = _fespace1.refFE().nbDofPerEdge();   // number of Dof per edge
    //UInt nDofPerFace1  = _fespace1.refFE().nbDofPerFace;   // number of Dof per face
    //UInt nDofPerElem1  = _fespace1.refFE().nbDofPerVolume; // number of Dof per Volume

    UInt nDofPerVert2  = _fespace2.refFE().nbDofPerVertex(); // number of Dof per vertex
    UInt nDofPerEdge2  = _fespace2.refFE().nbDofPerEdge();   // number of Dof per edge
    //UInt nDofPerFace2  = _fespace2.refFE().nbDofPerFace;   // number of Dof per face
    //UInt nDofPerElem2  = _fespace2.refFE().nbDofPerVolume; // number of Dof per Volume

    //UInt nElemV        = GeoShape::S_numVertices; // Number of element's vertices
    //UInt nElemE        = GeoShape::S_numEdges;    // Number of element's edges
    //UInt nElemF        = GeoShape::S_numFaces;    // Number of element's faces

    //UInt nBFacesVert   = GeoShape::GeoBShape::S_numVertices;
    //UInt nBFacesEdge   = GeoShape::GeoBShape::S_numEdges;
    //UInt nBFacesFace   = GeoShape::GeoBShape::S_numFaces;

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


    std::map<ID, ID> const& locDofMap = _dofInterface->localDofMap();
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
            if ((int)_fespace1.mesh()->pointList(iVert).marker() != M_data->fluidInterfaceFlag()) continue;

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
                    _vec2.setCoefficient( iDof, value );
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
            if ((int)_fespace1.mesh()->edgeList(iEdge).marker() != M_data->fluidInterfaceFlag()) continue;

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
                        _vec2.setCoefficient( iDof, value );
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
                        _vec2.setCoefficient( iDof1, value );
                        //
                        // ID of the 2nd node of the edge
                        ID node2 = _fespace1.mesh()->edgeList(iEdge).point(2).id();
                        iter = locDofMap.find(node2);
                        // ID of the 2nd dof of the edge in the solid numerotation
                        int iDof2 = icmp*numTotalDof2 + iter->second;
                        value = 0.5*_vec1( icmp*numTotalDof1 + iDofEdge ) + _vec2[iDof2];
                        //                                            std::cout << " transfering " << value << " from P2 " << iDofEdge << " to P1 " << iDof2 << " ( " << iter->second << " ) " << std::endl;
                        _vec2.setCoefficient( iDof2, value );
                        // now what to what boundary node ( in the solid numerotation ) should we send this value ?
                        //std::cout << "" << std::endl;
                        // Updating interpolated mesh velocity
                        //                                                     _vec2.setCoefficient( iDof2, value );
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
            if ((int)_fespace2.mesh()->edgeList(iEdge).marker() != M_data->fluidInterfaceFlag()) continue;
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
                    _vec2.setCoefficient(iDofEdge, value);
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
//                         _vec2.setCoefficient( iDof ,__sum);
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
//                 _vec2.setCoefficient( iDof, __sum);
//             }
//         }
}

} // Namespace LifeV

