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

/*!
    @file
    @brief Pure virtual operator class for FSI solvers

    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gilles Fourestey <fourestey@cscs.ch>
    @contributor Paolo Crosetto <paolo.crosetto@epfl.ch>
    @maintainer  Paolo Crosetto <paolo.crosetto@epfl.ch>

    @date 10-12-2010

    This is the base class for the FSI solvers in LifeV. It contains the methods to evaluate the residual and compute the
    Jacobian matrix, which make it suited for the generalized Newton algorithm implemented in NonlinearRichardson.hpp. The fluid
    and structure classes are members of this class and different formulations (e.g. Monolithic \cite CrosettoEtAl2009 , segregated
    Newton \cite FernandezMoubachir2005 , Dirichlet--Neumann \cite DeparisDiscacciati2006 , Robin Neumann \cite BadiaNobileVergara2008 ) are implemented in the derived classes.
 */

#include <lifev/core/LifeV.hpp>
#include <lifev/fsi/solver/FSIOperator.hpp>

namespace LifeV
{

// ===================================================
//  Constructors & Destructors
// ===================================================
FSIOperator::FSIOperator() :
    M_uFESpace                           ( ),
    M_pFESpace                           ( ),
    M_dFESpace                           ( ),
    M_mmFESpace                          ( ),
    M_fluidMesh                          ( ),
    M_solidMesh                          ( ),
    M_fluidLocalMesh                     ( ),
    M_solidLocalMesh                     ( ),
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
    M_fluidTimeAdvanceMethod             ( ),
    M_solidTimeAdvanceMethod             ( ),
    M_ALETimeAdvanceMethod               ( ),
    M_fluidTimeAdvance                   ( ),
    M_fluidMassTimeAdvance               ( ),
    M_solidTimeAdvance                   ( ),
    M_ALETimeAdvance                     ( ),
    M_dataFile                           ( ),
    M_meshDataFluid                      ( new MeshData() ),
    M_meshDataSolid                      ( new MeshData() ),
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
    M_rhs                                ( ),
    M_alphaF                             ( ), //vector_Type, for alphaf robin
    M_alphaFCoef                         ( 0 ),
    M_betaMean                           ( ),
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
    M_solidLeader                        ( ),
    M_aleOrder                           ( std::string ("P1") ),
    M_structureNonLinear                 (false)
{
}

FSIOperator::~FSIOperator()
{
}


// ===================================================
//  Virtual Public Methods
// ===================================================
void
FSIOperator::setDataFile ( const GetPot& dataFile )
{
    M_fluidMesh.reset ( new mesh_Type ( M_epetraComm ) );
    M_meshDataFluid->setup (dataFile, "fluid/space_discretization");
    readMesh (*M_fluidMesh, *M_meshDataFluid);

    M_solidMesh.reset (new mesh_Type ( M_epetraComm ) );
    M_meshDataSolid->setup (dataFile, "solid/space_discretization");
    readMesh (*M_solidMesh, *M_meshDataSolid);

    M_dataFile = dataFile;

    M_fluidTimeAdvanceMethod =  dataFile ( "fluid/time_discretization/method", "BDF");
    M_aleOrder               =  dataFile ( "fluid/space_discretization/ale_order", "P1");
    M_solidTimeAdvanceMethod =  dataFile ( "solid/time_discretization/method", "BDF");
    M_ALETimeAdvanceMethod   = dataFile ("mesh_motion/time_discretization/method", "BDF");
    M_structureNonLinear     = M_data->dataSolid()->solidType().compare ("linearVenantKirchhoff");
    this->setupTimeAdvance ( dataFile );
}

void
FSIOperator::setupFEspace()
{
    Displayer disp (M_epetraComm);
    disp.leaderPrint ("FSI-  Setting ReferenceFE and QuadratureRule ...           \n");

    std::string uOrder = M_data->dataFluid()->uOrder();
    std::string pOrder = M_data->dataFluid()->pOrder();
    std::string dOrder = M_data->dataSolid()->order();

    ASSERT ( ! ( (!uOrder.compare ("P2") && !dOrder.compare ("P1") ) ), "You are using P2 FE for the fluid velocity and P1 FE for the structure: when the coupling is enforced only in the P1 nodes, the energy is not conserved across the interface." );
    ASSERT ( ! ( (!uOrder.compare ("P1") && !dOrder.compare ("P2") ) ), "You are using P1 FE for the fluid velocity and P2 FE for the structure: when the coupling is enforced only in the P1 nodes, the energy is not conserved across the interface." );
    ASSERT ( ! ( (!uOrder.compare ("P1Bubble") && !dOrder.compare ("P2") ) ), "You are using P2 FE for the fluid velocity and P1 FE for the structure: when the coupling is enforced only in the P1 nodes, the energy is not conserved across the interface." );

    if ( M_aleOrder.compare ( M_meshDataFluid->mOrder() ) )
    {
        disp.leaderPrint ("FSI-  WARNING! The mesh order is different\n");
        disp.leaderPrint ("      => The nodes of the mesh will not entirely follow the computed displacement.\n");
    }

    const ReferenceFE*    refFE_vel (0);
    const QuadratureRule* qR_vel (0);
    const QuadratureRule* bdQr_vel (0);

    const ReferenceFE*    refFE_press (0);
    const QuadratureRule* qR_press (0);
    const QuadratureRule* bdQr_press (0);

    const ReferenceFE*    refFE_struct (0);
    const QuadratureRule* qR_struct (0);
    const QuadratureRule* bdQr_struct (0);

    const ReferenceFE*    refFE_mesh (0);
    const QuadratureRule* qR_mesh (0);
    const QuadratureRule* bdQr_mesh (0);

    if ( uOrder.compare ("P2") == 0 )
    {
        refFE_vel = &feTetraP2;
        qR_vel    = &quadRuleTetra15pt; // DoE 5
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( uOrder.compare ("P1") == 0 )
    {
        refFE_vel = &feTetraP1;
        qR_vel    = &quadRuleTetra4pt;  // DoE 2
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( uOrder.compare ("P1Bubble") == 0 )
    {
        refFE_vel = &feTetraP1bubble;
        qR_vel    = &quadRuleTetra64pt;  // DoE 2
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else
    {
        ERROR_MSG (uOrder + " velocity FE not implemented yet.");
    }

    if ( pOrder.compare ("P2") == 0 )
    {
        refFE_press = &feTetraP2;
        qR_press    = qR_vel; // DoE 5
        bdQr_press  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( pOrder.compare ("P1") == 0 )
    {
        refFE_press = &feTetraP1;
        qR_press    = qR_vel;  // DoE 2
        bdQr_press  = &quadRuleTria3pt;   // DoE 2
    }
    else
    {
        ERROR_MSG (pOrder + " pressure FE not implemented yet.");
    }

    if ( dOrder.compare ("P2") == 0 )
    {
        refFE_struct = &feTetraP2;
        qR_struct    = &quadRuleTetra15pt; // DoE 5
        bdQr_struct  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( dOrder.compare ("P1") == 0 )
    {
        refFE_struct = &feTetraP1;
        qR_struct    = &quadRuleTetra4pt;  // DoE 2
        bdQr_struct  = &quadRuleTria3pt;   // DoE 2
    }
    else
    {
        ERROR_MSG (dOrder + " structure FE not implemented yet.");
    }

    if ( M_aleOrder.compare ("P2") == 0 )
    {
        refFE_mesh = &feTetraP2;
        qR_mesh    = &quadRuleTetra15pt; // DoE 5
        bdQr_mesh  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( M_aleOrder.compare ("P1") == 0 )
    {
        refFE_mesh = &feTetraP1;
        qR_mesh    = &quadRuleTetra4pt;  // DoE 2
        bdQr_mesh  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( M_aleOrder.compare ("P1Bubble") == 0 )
    {
        refFE_mesh = &feTetraP1bubble;
        qR_mesh    = &quadRuleTetra64pt; // DoE 2
        bdQr_mesh  = &quadRuleTria3pt;   // DoE 2
    }
    else
    {
        ERROR_MSG (M_aleOrder + " mesh FE not implemented yet.");
    }

    disp.leaderPrint ("done\n");

    disp.leaderPrint ("FSI-  Building fluid FESpace ...               \n");
    if (this->isFluid() )
    {

        M_mmFESpace.reset (new FESpace<mesh_Type, MapEpetra> (M_fluidLocalMesh,
                                                              //dOrder,
                                                              *refFE_mesh,
                                                              *qR_mesh,
                                                              *bdQr_mesh,
                                                              3,
                                                              M_epetraComm) );

        M_uFESpace.reset ( new FESpace<mesh_Type, MapEpetra> (M_fluidLocalMesh,
                                                              //uOrder,
                                                              *refFE_vel,
                                                              *qR_vel,
                                                              *bdQr_vel,
                                                              3,
                                                              M_epetraComm) );

        M_pFESpace.reset ( new FESpace<mesh_Type, MapEpetra> (M_fluidLocalMesh,
                                                              //pOrder,
                                                              *refFE_press,
                                                              *qR_press,
                                                              *bdQr_press,
                                                              1,
                                                              M_epetraComm) );
    }
    else
    {
        M_mmFESpace.reset (new FESpace<mesh_Type, MapEpetra> (M_fluidMesh,
                                                              //dOrder,
                                                              *refFE_mesh,
                                                              *qR_mesh,
                                                              *bdQr_mesh,
                                                              3,
                                                              M_epetraComm) );

        M_uFESpace.reset ( new FESpace<mesh_Type, MapEpetra> (M_fluidMesh,
                                                              //uOrder,
                                                              *refFE_vel,
                                                              *qR_vel,
                                                              *bdQr_vel,
                                                              3,
                                                              M_epetraComm) );

        M_pFESpace.reset ( new FESpace<mesh_Type, MapEpetra> (M_fluidMesh,
                                                              //pOrder,
                                                              *refFE_press,
                                                              *qR_press,
                                                              *bdQr_press,
                                                              1,
                                                              M_epetraComm) );
    }
    M_epetraWorldComm->Barrier();

    disp.leaderPrint ("FSI-  Building solid FESpace ...               \n");
    if (this->isSolid() )
    {
        M_dFESpace.reset ( new FESpace<mesh_Type, MapEpetra> ( M_solidLocalMesh,
                                                               dOrder,
                                                               //*refFE_struct,
                                                               //*qR_struct,
                                                               //*bdQr_struct,
                                                               3,
                                                               M_epetraComm) );
    }
    else
    {
        M_dFESpace.reset (new FESpace<mesh_Type, MapEpetra> (M_solidMesh,
                                                             //dOrder,
                                                             *refFE_struct,
                                                             *qR_struct,
                                                             *bdQr_struct,
                                                             3,
                                                             M_epetraComm) );
    }
    M_epetraWorldComm->Barrier();
}


void
FSIOperator::partitionMeshes()
{
    if (this->isFluid() )
    {
        MeshPartitioner< mesh_Type > fluidPartitioner (M_fluidMesh, M_epetraComm);
        M_fluidLocalMesh = fluidPartitioner.meshPartition();
    }
    if (this->isSolid() )
    {
        MeshPartitioner< mesh_Type > solidPartitioner (M_solidMesh, M_epetraComm);
        M_solidLocalMesh = solidPartitioner.meshPartition();
    }
}

#ifdef HAVE_HDF5
void
FSIOperator::partitionMeshes ( meshFilter_Type& fluidMeshFilter, meshFilter_Type& solidMeshFilter )
{
    M_fluidMesh = fluidMeshFilter.getMeshPartition();
    M_solidMesh = solidMeshFilter.getMeshPartition();
}
#endif


void
FSIOperator::setupDOF ( void )
{
    Displayer disp (M_epetraWorldComm);
    disp.leaderPrint ("FSI: setting DOF ... " );
    DOF uDof (*M_fluidMesh, M_uFESpace->refFE() ); // velocity dof related to unpartitioned mesh
    DOF dDof (*M_solidMesh, M_dFESpace->refFE() ); // velocity dof related to unpartitioned mesh

    M_dofFluidToStructure                .reset ( new DOFInterface3Dto3D );
    M_dofStructureToFluid                .reset ( new DOFInterface3Dto3D );
    M_dofStructureToSolid                .reset ( new DOFInterface3Dto3D );
    M_dofStructureToHarmonicExtension    .reset ( new DOFInterface3Dto3D );
    M_dofHarmonicExtensionToFluid        .reset ( new DOFInterface3Dto3D );
    M_dofFluid                           .reset ( new DOFInterface3Dto2D );
    M_dofSolid                           .reset ( new DOFInterface3Dto2D );
    M_dofFluidInv                        .reset ( new DOFInterface3Dto2D );
    M_dofSolidInv                        .reset ( new DOFInterface3Dto2D );
    M_bcvFluidInterfaceDisp              .reset ( new  BCVectorInterface );
    M_bcvFluidLoadToStructure            .reset ( new  BCVectorInterface );
    M_bcvSolidLoadToStructure            .reset ( new  BCVectorInterface );
    M_bcvStructureToFluid                .reset ( new  BCVectorInterface );
    M_bcvStructureDispToFluid            .reset ( new  BCVectorInterface );
    M_bcvStructureDispToSolid            .reset ( new  BCVectorInterface );
    M_bcvStructureDispToHarmonicExtension.reset ( new  BCVectorInterface );
    M_bcvHarmonicExtensionVelToFluid     .reset ( new  BCVectorInterface );
    M_bcvDerHarmonicExtensionVelToFluid  .reset ( new  BCVectorInterface );
    M_bcvDerFluidLoadToStructure         .reset ( new  BCVectorInterface );
    M_bcvDerFluidLoadToFluid             .reset ( new  BCVectorInterface );
    M_bcvDerStructureDispToSolid         .reset ( new  BCVectorInterface );

    M_dofFluidToStructure->setup (   M_dFESpace->refFE(), M_dFESpace->dof(),
                                     M_uFESpace->refFE(), uDof );
    M_dofFluidToStructure->update ( *M_dFESpace->mesh(),  M_data->structureInterfaceFlag(),
                                    *M_fluidMesh,         M_data->fluidInterfaceFlag(),
                                    M_data->interfaceTolerance(),
                                    M_data->structureInterfaceVertexFlag() );

    //here the solid mesh must be non partitioned in the monolithic case
    M_dofStructureToHarmonicExtension->setup (   M_uFESpace->refFE(), M_uFESpace->dof(),
                                                 M_dFESpace->refFE(), dDof );
    M_dofStructureToHarmonicExtension->update ( *M_uFESpace->mesh(),  M_data->fluidInterfaceFlag(),
                                                *M_solidMesh, M_data->structureInterfaceFlag(),
                                                M_data->interfaceTolerance(),
                                                M_data->fluidInterfaceVertexFlag() );

    M_dofStructureToSolid->setup (   M_dFESpace->refFE(), M_dFESpace->dof(),
                                     M_dFESpace->refFE(), M_dFESpace->dof() );
    M_dofStructureToSolid->update ( *M_dFESpace->mesh(),  M_data->structureInterfaceFlag(),
                                    *M_dFESpace->mesh(),  M_data->structureInterfaceFlag(),
                                    M_data->interfaceTolerance(),
                                    M_data->structureInterfaceVertexFlag() );

    M_dofStructureToFluid->setup (   M_uFESpace->refFE(), M_uFESpace->dof(), //modifica matteo FSI
                                     M_dFESpace->refFE(), dDof );
    M_dofStructureToFluid->update ( *M_uFESpace->mesh(), M_data->fluidInterfaceFlag(),
                                    *M_solidMesh,        M_data->structureInterfaceFlag(),
                                    M_data->interfaceTolerance(),
                                    M_data->fluidInterfaceVertexFlag() );

    M_dofHarmonicExtensionToFluid->setup (   M_uFESpace->refFE(),  M_uFESpace->dof(),
                                             M_uFESpace->refFE(),  M_uFESpace->dof() );
    M_dofHarmonicExtensionToFluid->update ( *M_uFESpace->mesh(),  M_data->fluidInterfaceFlag(),
                                            *M_uFESpace->mesh(),  M_data->fluidInterfaceFlag(),
                                            M_data->interfaceTolerance(),
                                            M_data->fluidInterfaceVertexFlag() );

    M_epetraWorldComm->Barrier();
    disp.leaderPrint ("done\n");

    createInterfaceMaps ( M_dofStructureToHarmonicExtension->localDofMap() );
}

void
FSIOperator::setupFluidSolid ( void )
{
    setupFluidSolid (imposedFluxes() );
}

void
FSIOperator::setupFluidSolid ( UInt const fluxes )
{
    M_meshMotion.reset ( new meshMotion_Type ( *M_mmFESpace, M_epetraComm ) );

    if ( this->isFluid() )
    {
        M_fluid.reset ( new fluid_Type ( M_data->dataFluid(), *M_uFESpace, *M_pFESpace, M_epetraComm, fluxes ) );

        //Vector initialization
        M_rhs.reset ( new vector_Type ( M_fluid->getMap() ) );
    }

    M_solid.reset ( new solid_Type( ) );
    M_solid->setup ( M_data->dataSolid(), M_dFESpace, M_epetraComm );

    M_epetraWorldComm->Barrier();
}



void
FSIOperator::setupSystem ( void )
{
    if ( this->isFluid() )
    {
        //Data
        M_fluid->setUp ( M_dataFile );
        M_meshMotion->setUp ( M_dataFile );
        //         if (M_linearFluid)
        //             M_fluidLin->setUp(dataFile);
    }

    if ( this->isSolid() )
    {
        M_solid->setDataFromGetPot ( M_dataFile );
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

    if (this->isSolid() )
    {
        M_data->dataSolid()->showMe();
        //initialize xi0 for timaAdvance method for solid
        double  xi = M_solidTimeAdvance->coefficientSecondDerivative ( 0 ) / ( M_data->dataSolid()->dataTime()->timeStep() * M_data->dataSolid()->dataTime()->timeStep() );
        M_solid->buildSystem (xi);
        M_solid->massMatrix()->globalAssemble();
    }

}

void
FSIOperator::updateSystem()
{
    if ( this->isFluid() )
    {
        M_ALETimeAdvance->updateRHSContribution (M_data->dataFluid()->dataTime()->timeStep() );
        M_fluidTimeAdvance->updateRHSContribution (M_data->dataFluid()->dataTime()->timeStep() );

        if (M_data->dataFluid()->conservativeFormulation() )
        {
            M_fluidMassTimeAdvance->updateRHSContribution (M_data->dataFluid()->dataTime()->timeStep() );
        }

        transferMeshMotionOnFluid (M_meshMotion->disp(), *this->M_dispFluidMeshOld);

    }

    if ( this->isSolid() )
    {
        M_solid->updateSystem();
        M_solidTimeAdvance->updateRHSContribution ( M_data->dataSolid()->dataTime()->timeStep() );
        vector_Type rhsW (M_dFESpace->map(), Repeated);
        rhsW = (*M_solid->massMatrix() *  (M_solidTimeAdvance->rhsContributionSecondDerivative() ) * M_data->dataSolid()->dataTime()->timeStep() * M_data->dataSolid()->dataTime()->timeStep() / M_solidTimeAdvance->coefficientSecondDerivative ( 0 ) );
        M_solid->setRightHandSide ( rhsW ); //for the solid rhs;
    }

    couplingVariableExtrap( );

}



void FSIOperator::couplingVariableExtrap( )
{
    *M_lambda      = lambdaSolid();

    if (!M_lambdaDot.get() )
    {
        M_lambdaDot.reset   ( new vector_Type (*M_fluidInterfaceMap, Unique) );
        *M_lambda     += M_data->dataFluid()->dataTime()->timeStep() * lambdaDotSolid();
    }
    else
    {
        if ( this->isSolid() )
        {
            vector_Type solidDisp (M_solid->displacement() );
            vector_Type solidVel (M_solid->displacement() );
            M_solidTimeAdvance->extrapolation (solidDisp);
            M_solidTimeAdvance->extrapolationFirstDerivative (solidVel);
            transferSolidOnInterface (solidDisp, *M_lambda);
            transferSolidOnInterface (solidVel, *M_lambdaDot);
        }
    }
    displayer().leaderPrint ("FSI-  norm( disp  ) init =                     ", M_lambda->normInf(), "\n" );
    displayer().leaderPrint ("FSI-  norm( velo )  init =                     ", M_lambdaDot->normInf(), "\n");
}

void
FSIOperator::updateSolution ( const vector_Type& /* solution */ )
{
    if ( this->isFluid() )
    {
        M_ALETimeAdvance->shiftRight ( this->M_meshMotion->disp() );
        M_fluidTimeAdvance->shiftRight ( *M_fluid->solution() );
        if (M_data->dataFluid()->conservativeFormulation() )
        {
            M_fluidTimeAdvance->shiftRight ( M_fluid->matrixMass() * (*M_fluid->solution() ) );
        }
    }
    if ( this->isSolid() )
    {
        M_solidTimeAdvance->shiftRight ( M_solid->displacement() );
    }
}

UInt
FSIOperator::imposedFluxes ( void )
{
    if ( this->isFluid() )
    {
        std::vector<bcName_Type> fluxVector = M_BCh_u->findAllBCWithType ( Flux );
        UInt numLM = static_cast<UInt> ( fluxVector.size() );

        UInt offset = M_uFESpace->map().map (Unique)->NumGlobalElements()
                      + M_pFESpace->map().map (Unique)->NumGlobalElements();

        for ( UInt i = 0; i < numLM; ++i )
        {
            M_BCh_u->setOffset ( fluxVector[i], offset + i );
        }
        return numLM;
    }
    else
    {
        return 0;
    }
}

void
FSIOperator::initialize ( fluidPtr_Type::value_type::function_Type const& u0,
                          fluidPtr_Type::value_type::function_Type const& p0,
                          solidPtr_Type::value_type::function const& d0,
                          solidPtr_Type::value_type::function const& w0,
                          fluidPtr_Type::value_type::function_Type const& /*df0*/ )
{
    debugStream ( 6220 ) << "FSI:: solid init \n";
    if (this->isSolid() )
    {
        solid().initialize (d0, w0, w0);
    }
    debugStream ( 6220 ) << "FSI:: fluid init \n";
    if (this->isFluid() )
    {
        fluid().initialize (u0, p0);
    }
}

void
FSIOperator::setupTimeAdvance ( const dataFile_Type& dataFile )
{

    if (this->isFluid() )
    {
        M_fluidTimeAdvance.reset ( TimeAdvanceFactory::instance().createObject ( M_fluidTimeAdvanceMethod ) );
        if (M_data->dataFluid()->conservativeFormulation() )
        {
            M_fluidMassTimeAdvance.reset ( TimeAdvanceFactory::instance().createObject ( M_fluidTimeAdvanceMethod ) );
        }


        if (M_fluidTimeAdvanceMethod == "Newmark")
        {
            M_fluidTimeAdvance->setup (M_data->dataFluid()->dataTimeAdvance()->coefficientsNewmark() , 1);
            if (M_data->dataFluid()->conservativeFormulation() )
            {
                M_fluidMassTimeAdvance->setup (M_data->dataFluid()->dataTimeAdvance()->coefficientsNewmark() , 1);
            }
        }
        if (M_fluidTimeAdvanceMethod == "BDF")
        {
            M_fluidTimeAdvance->setup ( M_data->dataFluid()->dataTimeAdvance()->orderBDF(), 1 );
            if (M_data->dataFluid()->conservativeFormulation() )
            {
                M_fluidMassTimeAdvance->setup ( M_data->dataFluid()->dataTimeAdvance()->orderBDF(), 1 );
            }
        }
        M_ALETimeAdvance.reset ( TimeAdvanceFactory::instance().createObject ( M_ALETimeAdvanceMethod ) );

        if (M_ALETimeAdvanceMethod == "Newmark")
        {
            M_ALETimeAdvance->setup ( M_data->timeAdvanceDataALE()->coefficientsNewmark() , 2);
        }

        if (M_ALETimeAdvanceMethod == "BDF")
        {
            M_ALETimeAdvance->setup ( M_data->timeAdvanceDataALE()->orderBDF(), 1 );
        }

        M_fluidTimeAdvance->setTimeStep ( M_data->dataFluid()->dataTime()->timeStep() );
        if (M_data->dataFluid()->conservativeFormulation() )
        {
            M_fluidMassTimeAdvance->setTimeStep ( M_data->dataFluid()->dataTime()->timeStep() );
        }

        M_ALETimeAdvance->setTimeStep ( M_data->dataFluid()->dataTime()->timeStep() );

        if (this->isLeader() )
        {
            M_fluidTimeAdvance->showMe();
            //M_fluidMassTimeAdvance->showMe();
            M_ALETimeAdvance->showMe();
        }
    }
    if ( this->isSolid() )
    {
        M_solidTimeAdvance.reset ( TimeAdvanceFactory::instance().createObject ( M_solidTimeAdvanceMethod ) );

        std::vector<Real> parameters (2);
        parameters[0]  = dataFile ("solid/time_discretization/theta", 0.25);
        parameters[1]  = dataFile ("solid/time_discretization/zeta", 0.5);
        UInt order = dataFile ("solid/time_discretization/BDF_order", 1);
        // Real rhoInfty = dataFile("solid/time_discretization/rhoInf", 1.);
        std::string type    = dataFile ("mesh_motion/time_discretization/typeOfGeneralizedAlpha", "HHT");


        if (M_solidTimeAdvanceMethod == "Newmark")
        {
            M_solidTimeAdvance->setup ( parameters, 2 );
        }

        if (M_solidTimeAdvanceMethod == "BDF")
        {
            M_solidTimeAdvance->setup ( order , 2);
        }

        M_solidTimeAdvance->setTimeStep ( M_data->dataSolid()->dataTime()->timeStep() );
        if (this->isLeader() )
        {
            M_solidTimeAdvance->showMe();
        }
    }

}
// ===================================================
//  Public Methods
// ===================================================
void FSIOperator::initializeTimeAdvance ( const std::vector<vectorPtr_Type>& initialFluidVel, const std::vector<vectorPtr_Type>& initialSolidDisp, const std::vector<vectorPtr_Type>&  initialFluidDisp)
{
    displayer().leaderPrint ("initializeTimeAdvance start");

    if ( this->isFluid() )
    {
        ASSERT ( initialFluidVel.size() == M_fluidTimeAdvance->size(), "The number of vectors for initializing the time scheme for the fluid velocity is not consistent with the discretization chosen" );

        if (M_data->dataFluid()->conservativeFormulation() )
        {
            ASSERT ( initialFluidVel.size() == M_fluidMassTimeAdvance->size(), "The number of vectors for initializing the time scheme for the fluid velocity is not consistent with the discretization chosen" );
        }

        ASSERT (initialFluidDisp.size() == M_ALETimeAdvance->size() , "The number of vectors for initializing the time discretization for the ALE map is not consistent with the discretization chosen");
        this->M_fluidTimeAdvance->setInitialCondition (initialFluidVel);

        if (M_data->dataFluid()->conservativeFormulation() )
        {
            this->M_fluidMassTimeAdvance->setInitialCondition (initialFluidVel);
        }


        this->M_ALETimeAdvance->setInitialCondition (initialFluidDisp);
    }
    if ( this->isSolid() )
    {
        ASSERT (initialSolidDisp.size() == M_solidTimeAdvance->size(), "The number of vectors for initializing the time scheme for the structure displacement is not consistent with the discretization chosen" );
        this->M_solidTimeAdvance->setInitialCondition (initialSolidDisp);
    }
    displayer().leaderPrint ("initializeTimeAdvance end");
}

void
FSIOperator::initializeFluid ( const vector_Type& velAndPressure,
                               const vector_Type& displacement )
{
    this->fluid().initialize ( velAndPressure );
    this->moveMesh ( displacement);
}

void
FSIOperator::initializeSolid ( vectorPtr_Type displacement,
                               vectorPtr_Type velocity )
{
    this->solid().initialize ( displacement, velocity);
}

void
FSIOperator::moveMesh ( const vector_Type& dep )
{
    displayer().leaderPrint ("FSI-  Moving the mesh ...                      ");
    M_fluidLocalMesh->meshTransformer().moveMesh (dep,  this->M_mmFESpace->dof().numTotalDof() );
    displayer().leaderPrint ( "done\n" );
    M_fluid->setRecomputeMatrix ( true );
}

void FSIOperator::createInterfaceMaps ( std::map<ID, ID> const& locDofMap )
{
    Displayer disp (M_epetraWorldComm);
    disp.leaderPrint ("FSI-  Building fluid variables ...             ");

    // now we build the sigma and lambda variables on each proc


    std::vector<int> dofInterfaceFluid;

    //is the interface map between HE (first) and solid (second)
    if ( this->isFluid() )
    {
        //std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
        dofInterfaceFluid.reserve ( locDofMap.size() );

        for (UInt dim = 0; dim < nDimensions; ++dim)
            for ( iterator_Type i = locDofMap.begin(); i != locDofMap.end(); ++i )
            {
                dofInterfaceFluid.push_back (i->second + dim * M_dFESpace->dof().numTotalDof() );    // in solid numerotation
            }
    }

    int* pointerToDofs (0);
    if (dofInterfaceFluid.size() > 0)
    {
        pointerToDofs = &dofInterfaceFluid[0];
    }

    M_fluidInterfaceMap.reset ( new MapEpetra ( -1,
                                                static_cast<int> (dofInterfaceFluid.size() ),
                                                pointerToDofs,
                                                M_epetraWorldComm ) );
    disp.leaderPrint ("done\n");
    M_epetraWorldComm->Barrier();

    disp.leaderPrint ("FSI-  Building solid variables ...             ");

    std::vector<int> dofInterfaceSolid;

    if (this->isSolid() )
    {
        //std::map<ID, ID> const& locDofMap = M_dofFluidToStructure->locDofMap();
        dofInterfaceSolid.reserve ( locDofMap.size() );
        //std::cout << "solid" << std::endl;
        for (UInt dim = 0; dim < nDimensions; ++dim)
            for ( iterator_Type i = locDofMap.begin(); i != locDofMap.end(); ++i )
            {
                dofInterfaceSolid.push_back (i->second/*Simone: first*/ + dim * M_dFESpace->dof().numTotalDof() ); // in solid numerotation
            }
    }


    pointerToDofs = 0;
    if (dofInterfaceSolid.size() > 0)
    {
        pointerToDofs = &dofInterfaceSolid[0];
    }

    M_solidInterfaceMap.reset ( new MapEpetra ( -1,
                                                static_cast<int> (dofInterfaceSolid.size() ),
                                                pointerToDofs,
                                                M_epetraWorldComm ) );

    M_epetraWorldComm->Barrier();
    disp.leaderPrint ("done\n");

    disp.leaderPrint ("FSI-  Variables initialization ...             \n");

    variablesInit ( M_data->dataSolid()->order() );

    M_epetraWorldComm->Barrier();
}

void
FSIOperator::transferFluidOnInterface (const vector_Type& _vec1, vector_Type& _vec2)
{
    // e.g.: vec1=M_fluid->residual(), vec2=M_sigmaFluid
    //     _vec2 = ZeroVector(_vec2.size());

    if (_vec1.mapType() == Unique)
    {
        vector_Type const  vec1Repeated (_vec1, Repeated);
        transferFluidOnInterface (vec1Repeated, _vec2);
        return;
    }

    if (_vec2.mapType() == Repeated)
    {
        vector_Type  vec2Unique (_vec2, Unique);
        transferFluidOnInterface (_vec1, vec2Unique);
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
            _vec2.setCoefficient ( it->second + dim * numTotalDofSolid,
                                   _vec1[it->first + dim * numTotalDofFluid] );
        }
}

//works in serial but no yet in parallel
void
FSIOperator::transferSolidOnFluid (const vector_Type& _vec1, vector_Type& _vec2)
{
    //    e.g.: vec2=M_fluid->residual(), vec1=M_sigmaFluid
    //     _vec2 = ZeroVector(_vec2.size());

    if (_vec1.mapType() == Unique)
    {
        vector_Type const  vec1Repeated (_vec1, Repeated);
        transferSolidOnInterface (vec1Repeated, _vec2);
        return;
    }

    if (_vec2.mapType() == Repeated)
    {
        vector_Type  vec2Unique (_vec2, Unique);
        transferSolidOnInterface (_vec1, vec2Unique);
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
            _vec2.setCoefficient ( it->second + dim * numTotalDofFluid,
                                   _vec1[it->first + dim * numTotalDofSolid]);
        }
}

void
FSIOperator::transferSolidOnInterface (const vector_Type& _vec1, vector_Type& _vec2)
{
    /* e.g.:
       vec1 (Unique)       vec2 (Repeated)  (On different MapEpetras) (Changing now to Unique on both)
       M_solid->disp()     M_lambdaSolid
       M_solid->vel()      M_lambdaDotSolid
       M_solid->residual() M_sigmaSolid
    */

    if (_vec1.mapType() == Unique)
    {
        vector_Type const  vec1Repeated (_vec1, Repeated);
        transferSolidOnInterface (vec1Repeated, _vec2);
        return;
    }

    if (_vec2.mapType() == Repeated)
    {
        vector_Type  vec2Unique (_vec2, Unique);
        transferSolidOnInterface (_vec1, vec2Unique);
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
            _vec2.setCoefficient ( it->second + dim * numTotalDofSolid,
                                   _vec1[it->first + dim * numTotalDofSolid] );
        }
}

void
FSIOperator::transferInterfaceOnSolid (const vector_Type& _vec1, vector_Type& _vec2)
{
    /* e.g.:
       vec2                vec1
       M_solid->disp()     M_lambdaSolid
       M_solid->vel()      M_lambdaDotSolid
       M_solid->residual() M_sigmaSolid
    */

    if (_vec1.mapType() == Unique)
    {
        vector_Type const  vec1Repeated (_vec1, Repeated);
        transferInterfaceOnSolid (vec1Repeated, _vec2);
        return;
    }

    if (_vec2.mapType() == Repeated)
    {
        vector_Type  vec2Unique (_vec2, Unique);
        transferInterfaceOnSolid (_vec1, vec2Unique);
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
            _vec2.setCoefficient ( it->second + dim * numTotalDofSolid,
                                   _vec1[it->second + dim * numTotalDofSolid] );
        }
}

void
FSIOperator::bcManageVectorRHS ( const fluidBchandlerPtr_Type& bcHandlerFluid, vector_Type& rhs )
{
    if ( !bcHandlerFluid->bcUpdateDone() || M_fluid->recomputeMatrix() )
    {
        bcHandlerFluid->bcUpdate ( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    }

    bcManageRhs ( rhs, *M_uFESpace->mesh(), M_uFESpace->dof(),  *bcHandlerFluid, M_uFESpace->feBd(), 0., 1. );
}

void
FSIOperator::bcManageVectorRHS ( const fluidBchandlerPtr_Type& bcHandlerFluid, const solidBchandlerPtr_Type& bcHandlerSolid, vector_Type& rhs )
{
    if ( !bcHandlerFluid->bcUpdateDone() || M_fluid->recomputeMatrix() )
    {
        bcHandlerFluid->bcUpdate ( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
    }

    bcManageRhs ( rhs, *M_uFESpace->mesh(), M_uFESpace->dof(),  *bcHandlerFluid, M_uFESpace->feBd(), 1., 0. );

    if ( !bcHandlerSolid->bcUpdateDone() )
    {
        bcHandlerSolid->bcUpdate ( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
    }

    bcManageRhs ( rhs, *M_dFESpace->mesh(), M_dFESpace->dof(),  *bcHandlerSolid, M_dFESpace->feBd(), 1., 0. );
}

void
FSIOperator::setAlphafCoef( )
{
    Real h = 0.1, R = 0.5;
    M_alphaFCoef  =  M_ALETimeAdvance->coefficientSecondDerivative ( 0 ) * (this->dataSolid()->rho() * h) / this->dataFluid()->dataTime()->timeStep();
    M_alphaFCoef += h * M_data->dataSolid()->young (1) * this->dataFluid()->dataTime()->timeStep() /
                    (pow (R, 2) * (1 - pow (M_data->dataSolid()->poisson (1), 2) ) );
    M_alphaFCoef /= M_ALETimeAdvance->coefficientFirstDerivative ( 0 );
}

void
FSIOperator::setStructureToFluidParameters()
{
    this->setAlphafCoef();
    this->setAlphaf();

    if (M_alphaF.get() == 0)
    {
        this->setAlphafCoef();
        M_bcvStructureToFluid->setRobinCoeff (M_alphaFCoef);
        M_bcvStructureToFluid->setBetaCoeff (M_alphaFCoef);
    }
    else
    {
        M_bcvStructureToFluid->setRobinCoeffVector (this->Alphaf() );
        M_bcvStructureToFluid->setBetaCoeffVector (this->Alphaf() );
    }
}

// ===================================================
//  Display Methods
// ===================================================
bool
FSIOperator::isLeader() const
{
    if ( isFluid() )
    {
        if ( M_fluid.get() == 0 )
        {
            return ( M_epetraComm->MyPID() == 0 );
        }

        return M_fluid->getDisplayer().isLeader();
    }

    if ( M_solid.get() == 0 )
    {
        return ( M_epetraComm->MyPID() == 0 );
    }

    return M_solid->displayer().isLeader();
}



Displayer const&
FSIOperator::displayer()
{
    if ( isFluid() &&  M_fluid.get() )
    {
        return M_fluid->getDisplayer();
    }

    if ( !isSolid() || !M_solid.get() )
    {
        std::cout << "ERROR: displayer not ready" << std::endl;
    }

    return M_solid->displayer();
}





// ===================================================
//  Get Functions
// ===================================================
/*
FSI::vector_Type
FSI::displacementOnInterface()
{

//     vector_Type dispOnInterface();
// //@     dispOnInterface = ZeroVector(dispOnInterface.size());

//  //    FOR_EACH_INTERFACE_DOF( dispOnInterface[IDsolid + jDim*totalDofSolid] =
// //                             M_solid->disp()[IDsolid + jDim*totalDofSolid]);

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
//  Set Functions
// ===================================================
void
FSIOperator::setComm ( const commPtr_Type& comm,
                       const commPtr_Type& worldComm )
{
    M_epetraComm       = comm;
    M_epetraWorldComm  = worldComm;
}


void
FSIOperator::setFluid ( const fluidPtr_Type& fluid, const meshMotionPtr_Type& meshmotion )
{
    M_fluid = fluid;
    M_meshMotion = meshmotion;
    M_isFluid = true;
}



void
FSIOperator::setSolid ( const solidPtr_Type& solid )
{
    M_solid = solid;
    M_isSolid = true;
}



void
FSIOperator::setFluidBC ( const fluidBchandlerPtr_Type& bc_fluid )
{
    if ( isFluid() )
    {
        M_BCh_u = bc_fluid;
    }
}



void
FSIOperator::setHarmonicExtensionBC ( const fluidBchandlerPtr_Type& bc_he )
{
    if ( isFluid() )
    {
        M_BCh_mesh = bc_he;
        //M_meshMotion->setBC(*M_BCh_mesh);
    }
}



void
FSIOperator::setSolidBC ( const solidBchandlerPtr_Type& bc_solid )
{
    if ( isSolid() )
    {
        M_BCh_d = bc_solid;
    }
}



void
FSIOperator::setLambdaFluid ( const vector_Type& lambda )
{
    if ( lambda.mapType() == Unique )
    {
        *M_lambdaFluid = lambda;
    }
    else // to be coded, I am not sure we need this functionality.
    {
        assert (false);    // if you get here, reformulate your problem in order to get a unique map as entry
    }

    *M_lambdaFluidRepeated = *M_lambdaFluid;
}



void
FSIOperator::setLambdaSolid ( const vector_Type& lambda )
{
    if ( lambda.mapType() == Unique )
    {
        *M_lambdaSolid = lambda;
    }
    else // to be coded, I am not sure we need this functionality.
    {
        assert (false);    // if you get here, reformulate your problem in order to get a unique map as entry
    }

    *M_lambdaSolidRepeated = *M_lambdaSolid;
}



void
FSIOperator::setLambdaSolidOld ( const vector_Type& lambda )
{
    if ( lambda.mapType() == Unique )
    {
        *M_lambdaSolidOld = lambda;
    }
    else // to be coded, I am not sure we need this functionality.
    {
        assert (false);    // if you get here, reformulate your problem in order to get a unique map as entry
    }
}



void
FSIOperator::setLambdaDotSolid ( const vector_Type& lambda )
{
    if ( lambda.mapType() == Unique )
    {
        *M_lambdaDotSolid = lambda;
    }
    else // to be coded, I am not sure we need this functionality.
    {
        assert (false);    // if you get here, reformulate your problem in order to get a unique map as entry
    }

    *M_lambdaDotSolidRepeated = *M_lambdaDotSolid;
}



void
FSIOperator::setSigmaSolid ( const vector_Type& sigma )
{
    if ( sigma.mapType() == Unique )
    {
        *M_sigmaSolid = sigma;
    }
    else // to be coded, I am not sure we need this functionality.
    {
        assert (false);    // if you get here, reformulate your problem in order to get a unique map as entry
    }

    *M_sigmaSolidRepeated = *M_sigmaSolid;
}


void
FSIOperator::setSigmaFluid ( const vector_Type& sigma )
{
    if ( sigma.mapType() == Unique )
    {
        *M_sigmaFluid = sigma;
    }
    else // to be coded, I am not sure we need this functionality.
    {
        assert (false);    // if you get here, reformulate your problem in order to get a unique map as entry
    }

    *M_sigmaFluidRepeated = *M_sigmaFluid;

}



void
FSIOperator::setMinusSigmaFluid ( const vector_Type& sigma )
{
    if ( sigma.mapType() == Unique )
    {
        *M_minusSigmaFluid  = sigma;
        *M_minusSigmaFluid *= -1;
    }
    else // to be coded, I am not sure we need this functionality.
    {
        assert (false);    // if you get here, reformulate your problem in order to get a unique map as entry
    }

    *M_minusSigmaFluidRepeated = *M_minusSigmaFluid;
}




void
FSIOperator::setAlphafbcf ( const bcFunction_Type& alphafbcf )
{
    vector_Type vec ( M_fluid->velocityFESpace().map() );
    M_fluid->velocityFESpace().interpolate ( static_cast<FESpace<mesh_Type, MapEpetra>::function_Type> ( alphafbcf ), vec, 0.0);
    *M_alphaF = vec ;
}



void
FSIOperator::setStructureDispToHarmonicExtension ( const vector_Type& disp, UInt type )
{
    M_bcvStructureDispToHarmonicExtension->setup ( disp,
                                                   M_dFESpace->dof().numTotalDof(),
                                                   M_dofStructureToHarmonicExtension,
                                                   type );
}



void
FSIOperator::setStructureToFluid ( const vector_Type& velo,  UInt type )
{
    M_bcvStructureToFluid->setup ( velo,
                                   M_uFESpace->dof().numTotalDof(),
                                   M_dofHarmonicExtensionToFluid,
                                   type );
}



void
FSIOperator::setStructureDispToFluid ( const vector_Type& disp,  UInt type )
{
    M_bcvStructureDispToFluid->setup ( disp,
                                       M_uFESpace->dof().numTotalDof(),
                                       M_dofStructureToFluid,
                                       type );
}



void
FSIOperator::setStructureDispToSolid ( const vector_Type& disp, UInt type )
{
    M_bcvStructureDispToSolid->setup ( disp,
                                       M_dFESpace->dof().numTotalDof(),
                                       M_dofStructureToSolid,
                                       type );
}



void
FSIOperator::setDerStructureDispToSolid ( const vector_Type& ddisp, UInt type )
{
    M_bcvDerStructureDispToSolid->setup ( ddisp,
                                          M_dFESpace->dof().numTotalDof(),
                                          M_dofStructureToSolid,
                                          type );
}



void
FSIOperator::setSolidLoadToStructure ( const vector_Type& load, UInt type )
{
    M_bcvSolidLoadToStructure->setup ( load,
                                       M_dFESpace->dof().numTotalDof(),
                                       M_dofStructureToFluid,
                                       type );
}



void
FSIOperator::setHarmonicExtensionVelToFluid ( const vector_Type& vel, UInt type )
{
    M_bcvHarmonicExtensionVelToFluid->setup ( vel,
                                              M_uFESpace->dof().numTotalDof(),
                                              M_dofHarmonicExtensionToFluid,
                                              type );
}



void
FSIOperator::setDerHarmonicExtensionVelToFluid ( const vector_Type& dvel, UInt type )
{
    M_bcvDerHarmonicExtensionVelToFluid->setup ( dvel,
                                                 M_uFESpace->dof().numTotalDof(),
                                                 M_dofHarmonicExtensionToFluid,
                                                 type );
}




void
FSIOperator::setFluidLoadToStructure ( const vector_Type& load, UInt type )
{
    M_bcvFluidLoadToStructure->setup ( load,
                                       M_dFESpace->dof().numTotalDof(),
                                       M_dofStructureToSolid,
                                       type );
}



void
FSIOperator::setDerFluidLoadToStructure ( const vector_Type& dload, UInt type )
{
    M_bcvDerFluidLoadToStructure->setup ( dload,
                                          M_dFESpace->dof().numTotalDof(),
                                          M_dofStructureToSolid,
                                          type );
}



void
FSIOperator::setDerFluidLoadToFluid ( const vector_Type& dload, UInt type )
{
    M_bcvDerFluidLoadToFluid->setup ( dload,
                                      M_uFESpace->dof().numTotalDof(),
                                      M_dofHarmonicExtensionToFluid,
                                      type );
}

void FSIOperator::setRobinOuterWall (function_Type const& dload, function_Type const& E)
{
    M_bcfRobinOuterWall.setFunctions_Robin (dload,
                                            E);
}


// ===================================================
//  Protected Methods
// ===================================================



void
FSIOperator::variablesInit ( const std::string& /*dOrder*/ )
//FSI::variablesInit(const ReferenceFE* refFE_struct,const LifeV::QuadratureRule*  bdQr_struct, const LifeV::QuadratureRule* qR_struct)
{
    M_lambdaFluid.reset        ( new vector_Type (*M_fluidInterfaceMap, Unique) );
    M_lambda.reset             ( new vector_Type (*M_solidInterfaceMap, Unique) );
    M_lambdaDot.reset             ( new vector_Type (*M_solidInterfaceMap, Unique) );
    M_lambdaFluidRepeated.reset ( new vector_Type (*M_fluidInterfaceMap, Repeated) );

    if ( this->isFluid() )
    {
        M_dispFluidMeshOld.reset ( new vector_Type (M_uFESpace->map(), Repeated) );
        M_veloFluidMesh.reset   ( new vector_Type (M_uFESpace->map(), Repeated) );
        M_alphaF.reset          ( new vector_Type (M_uFESpace->map(), Repeated) );

        if ( M_linearFluid )
        {
            M_derVeloFluidMesh.reset ( new vector_Type (this->M_uFESpace->map(), Repeated) );
        }
    }

    M_sigmaFluid.reset              ( new vector_Type (*M_fluidInterfaceMap, Unique) );
    M_sigmaFluidRepeated.reset      ( new vector_Type (*M_fluidInterfaceMap, Repeated) );
    M_minusSigmaFluid.reset         ( new vector_Type (*M_fluidInterfaceMap, Unique) );
    M_minusSigmaFluidRepeated.reset ( new vector_Type (*M_fluidInterfaceMap, Repeated) );

    M_lambdaSolid.reset   ( new vector_Type (*M_solidInterfaceMap, Unique) );
    M_lambdaSolidOld.reset ( new vector_Type (*M_solidInterfaceMap, Unique) );
    M_lambdaDotSolid.reset ( new vector_Type (*M_solidInterfaceMap, Unique) );
    M_sigmaSolid.reset    ( new vector_Type (*M_solidInterfaceMap, Unique) );

    M_lambdaSolidRepeated.reset   ( new vector_Type (*M_solidInterfaceMap, Repeated) );
    M_lambdaDotSolidRepeated.reset ( new vector_Type (*M_solidInterfaceMap, Repeated) );
    M_sigmaSolidRepeated.reset    ( new vector_Type (*M_solidInterfaceMap, Repeated) );
}


void
FSIOperator::transferMeshMotionOnFluid ( const vector_Type& _vec1, vector_Type& _vec2 )
{
    //transferMeshMotionOnFluid should handle the repetition of the interface nodes.
    if (_vec1.mapType() == Unique)
    {
        vector_Type const  vec1Repeated (_vec1, Repeated);
        transferMeshMotionOnFluid (vec1Repeated, _vec2);
        return;
    }

    if (_vec2.mapType() == Repeated)
    {
        vector_Type  vec2Unique (_vec2, Unique);
        transferMeshMotionOnFluid (_vec1, vec2Unique);
        _vec2 = vec2Unique;
        return;
    }

    _vec2 *= 0;

    interpolateVelocity (_vec1, _vec2);
    return;

}

void
FSIOperator::interpolateVelocity ( const vector_Type& _vec1, vector_Type& _vec2 )
{
    assert (_vec1.mapType() == Repeated);
    assert (_vec2.mapType() == Unique);

    typedef mesh_Type::elementShape_Type GeoShape; // Element shape

    UInt nDofpV = M_uFESpace->refFE().nbDofPerVertex(); // number of DOF per vertex
    UInt nDofpE = M_uFESpace->refFE().nbDofPerEdge();   // number of DOF per edge
    UInt nDofpF = M_uFESpace->refFE().nbDofPerFace();   // number of DOF per face
    UInt nDofpEl = M_uFESpace->refFE().nbDofPerVolume(); // number of DOF per Volume

    UInt nElemV = GeoShape::S_numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::S_numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::S_numFaces;    // Number of element's faces

    //    UInt nDofElem = M_uFESpace->refFE().nbDof; // Number of local dof per element of the M_uFESpace->mesh() (_mesh.getRefFE().nbDof)
    UInt nDofElemMesh = M_mmFESpace->refFE().nbDof();

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's DOF on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's DOF on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's DOF on a Element

    Real x, y, z;
    Vector wLoc ( nDofElemMesh * nDimensions );
    ID lDof;

    // Loop on elements of the mesh
    for ( ID iElem = 0; iElem < M_uFESpace->mesh()->numVolumes(); ++iElem )
    {
        UInt elemId = M_uFESpace->mesh()->volume ( iElem ).localId();
        //if (elemId != iElem)
        // std::cout << " elemId = " << elemId << " iElem = " << iElem << std::endl;

        // Updating the local mesh velocity in this mesh elment
        for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
            for ( ID idof = 0; idof < nDofElemMesh; ++idof )
            {
                wLoc ( icmp * nDofElemMesh + idof ) =
                    _vec1 ( icmp * M_mmFESpace->dof().numTotalDof()
                            + M_mmFESpace->dof().localToGlobalMap ( iElem, idof) );
            }
        // Vertex based Dof
        if ( nDofpV )
        {

            // loop on element vertices
            for ( ID iVe = 0; iVe < nElemV; ++iVe )
            {

                // Loop number of DOF per vertex
                for ( ID l = 0; l < nDofpV; ++l )
                {
                    lDof = iVe * nDofpV + l; // Local dof in this element

                    // Nodal coordinates
                    x = M_uFESpace->refFE().xi ( lDof );
                    y = M_uFESpace->refFE().eta ( lDof );
                    z = M_uFESpace->refFE().zeta ( lDof );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        Real __sum = 0;
                        for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local DOF on the element
                        {
                            __sum += wLoc ( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi ( idof, x, y, z );
                        }

                        // Updating interpolated mesh velocity
                        int iDof = icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobalMap ( iElem, lDof  );
                        _vec2.setCoefficient ( iDof , __sum);

                    }
                }
            }
        }

        // Edge based Dof
        if ( nDofpE )
        {

            // loop on element edges
            for ( ID iEd = 0; iEd < nElemE; ++iEd )
            {

                // Loop number of DOF per edge
                for ( ID l = 0; l < nDofpE; ++l )
                {
                    lDof = nDofElemV + iEd * nDofpE + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    x = M_uFESpace->refFE().xi ( lDof );
                    y = M_uFESpace->refFE().eta ( lDof );
                    z = M_uFESpace->refFE().zeta ( lDof );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        Real __sum = 0;
                        for ( ID idof = 0; idof < nDofElemMesh; ++idof )   // Loop on local DOF on the adjacent element
                        {
                            __sum += wLoc ( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi ( idof, x, y, z );    // Problem here with P2
                        }

                        // Updating interpolating vector
                        int iDof = icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobalMap ( iElem, lDof );
                        _vec2.setCoefficient ( iDof , __sum);

                    }
                }
            }
        }

        // Face based Dof
        if ( nDofpF )
        {

            // loop on element faces
            for ( ID iFa = 0; iFa < nElemF; ++iFa )
            {

                // Loop on number of DOF per face
                for ( ID l = 0; l < nDofpF; ++l )
                {

                    lDof = nDofElemE + nDofElemV + iFa  * nDofpF + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    x = M_uFESpace->refFE().xi ( lDof );
                    y = M_uFESpace->refFE().eta ( lDof );
                    z = M_uFESpace->refFE().zeta ( lDof );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        Real __sum = 0;
                        for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local DOF on the adjacent element
                        {
                            __sum += wLoc ( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi ( idof, x, y, z );    // Problem here with P2
                        }

                        // Updating interpolating vector
                        int iDof = icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobalMap ( iElem, lDof );
                        _vec2.setCoefficient ( iDof , __sum);
                    }
                }
            }
        }

        // Element based Dof
        // Loop on number of DOF per Element
        for ( ID l = 0; l < nDofpEl; ++l )
        {
            lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

            // Nodal coordinates
            x = M_uFESpace->refFE().xi ( lDof );
            y = M_uFESpace->refFE().eta ( lDof );
            z = M_uFESpace->refFE().zeta ( lDof );

            // Loop on data vector components
            for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
            {

                // Interpolating data at the nodal point
                Real __sum = 0;
                for ( ID idof = 0; idof < nDofElemMesh; ++idof )  // Loop on local DOF on the adjacent element
                {
                    __sum += wLoc ( icmp * nDofElemMesh + idof ) * M_mmFESpace->refFE().phi ( idof, x, y, z );
                }

                // Updating interpolating vector
                //                std::cout << M_uFESpace->dof().localToGlobal( elemId, lDof ) << " ";
                //                std::cout << icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobal( elemId, lDof ) << std::endl;

                int iDof = icmp * M_uFESpace->dof().numTotalDof() + M_uFESpace->dof().localToGlobalMap ( elemId, lDof );
                _vec2.setCoefficient ( iDof, __sum);
            }
        }
    }

}




// this will interpolate dofs values from fespace1 to fespace2
void
FSIOperator::interpolateInterfaceDofs ( const FESpace<mesh_Type, MapEpetra>& _fespace1,
                                        const vector_Type&                   _vec1,
                                        const FESpace<mesh_Type, MapEpetra>& _fespace2,
                                        vector_Type&                         _vec2,
                                        dofInterface3DPtr_Type&                _dofInterface)
{
    assert (_vec1.mapType() == Repeated);
    assert (_vec2.mapType() == Unique);

    typedef mesh_Type::elementShape_Type GeoShape; // Element shape


    UInt nDofPerVert1  = _fespace1.refFE().nbDofPerVertex(); // number of DOF per vertex
    UInt nDofPerEdge1  = _fespace1.refFE().nbDofPerEdge();   // number of DOF per edge
    //UInt nDofPerFace1  = _fespace1.refFE().nbDofPerFace;   // number of DOF per face
    //UInt nDofPerElem1  = _fespace1.refFE().nbDofPerVolume; // number of DOF per Volume

    UInt nDofPerVert2  = _fespace2.refFE().nbDofPerVertex(); // number of DOF per vertex
    UInt nDofPerEdge2  = _fespace2.refFE().nbDofPerEdge();   // number of DOF per edge
    //UInt nDofPerFace2  = _fespace2.refFE().nbDofPerFace;   // number of DOF per face
    //UInt nDofPerElem2  = _fespace2.refFE().nbDofPerVolume; // number of DOF per Volume

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

    //UInt nDofElemVert1 = nElemV * nDofPerVert1; // number of vertex's DOF on a Element
    //UInt nDofElemEdge1 = nElemE * nDofPerEdge1; // number of edge's DOF on a Element
    //UInt nDofElemFace1 = nElemF * nDofPerFace1; // number of face's DOF on a Element

    //UInt nDofElemVert2 = nElemV * nDofPerVert2; // number of vertex's DOF on a Element
    //UInt nDofElemEdge2 = nElemE * nDofPerEdge2; // number of edge's DOF on a Element
    //UInt nDofElemFace2 = nElemF * nDofPerFace2; // number of face's DOF on a Element

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
        for ( ID iVert = 0; iVert < _fespace1.mesh()->numVertices(); ++iVert )
        {
            if ( (int) _fespace1.mesh()->pointList (iVert).markerID() != M_data->fluidInterfaceFlag() )
            {
                continue;
            }

            ID nodeID = _fespace1.mesh()->pointList (iVert).id();
            // Loop number of DOF per vertex
            //TODO check this, loop does not depend on l
            for ( ID l = 0; l < nDofPerVert1; ++l )
            {
                // Loop on data vector components
                for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                {
                    // "Interpolating" data at the nodal point
                    iter = locDofMap.find (nodeID);
                    Real value = _vec1 ( icmp * numTotalDof1 + nodeID );
                    // now what to what boundary node ( in the solid numerotation ) should we send this value ?
                    //std::cout << "" << std::endl;
                    int iDof = icmp * numTotalDof2 + iter->second;
                    //                                    std::cout << " transfering " << value << " from P1 " << nodeID << " to P1 " << iDof  << " ( " << iter->second << " ) " << std::endl;
                    // Updating interpolated mesh velocity
                    _vec2.setCoefficient ( iDof, value );
                }
            }
        }
    }

    // Edge based Dof
    if ( nDofPerEdge1 )
    {
        // loop on boundary edges
        for ( ID iEdge = 0; iEdge < nBEdges1; ++iEdge )
        {
            if ( (int) _fespace1.mesh()->edgeList (iEdge).markerID() != M_data->fluidInterfaceFlag() )
            {
                continue;
            }

            // edge ID
            ID edgeID = _fespace1.mesh()->edgeList (iEdge).id();
            // dof position of the edge since unknowns are store using [ node | edges | faces | volumes ]
            int iDofEdge = numTotalVert1 + edgeID;

            if (nDofPerEdge2) // the second FE space has dofs on its edges
            {
                for ( ID l = 0; l < nDofPerEdge1; ++l )
                {
                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                    {
                        // ID of the dof in the solid numerotation
                        iter = locDofMap.find (iDofEdge);
                        int iDof = icmp * numTotalDof2 + iter->second;
                        // "Interpolating" data at the nodal point
                        Real value = _vec1 ( icmp * numTotalDof1 + iDofEdge );
                        // now what to what boundary node ( in the solid numerotation ) should we send this value ?
                        //                                            std::cout << " transfering " << value << " from P2 " << edgeID << " to P2 " << iDof  << " ( " << iter->second << " ) " << std::endl;
                        // Updating interpolated mesh velocity
                        _vec2.setCoefficient ( iDof, value );
                    }
                }
            }
            else
            {
                for ( ID l = 0; l < nDofPerEdge1; ++l )
                {
                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                    {
                        //
                        // ID of the 1st node of the edge
                        ID node1 = _fespace1.mesh()->edgeList (iEdge).point (0).id();
                        iter = locDofMap.find (node1);
                        // ID of the 1st dof of the edge in the solid numerotation
                        int iDof1 = icmp * numTotalDof2 + iter->second;
                        value = 0.5 * _vec1 ( icmp * numTotalDof1 + iDofEdge  ) + _vec2[iDof1];
                        //                                            std::cout << " transfering " << value << " from P2 " << iDofEdge << " to P1 " << iDof1 << " ( " << iter->second << " ) " << std::endl;
                        _vec2.setCoefficient ( iDof1, value );
                        //
                        // ID of the 2nd node of the edge
                        ID node2 = _fespace1.mesh()->edgeList (iEdge).point (1).id();
                        iter = locDofMap.find (node2);
                        // ID of the 2nd dof of the edge in the solid numerotation
                        int iDof2 = icmp * numTotalDof2 + iter->second;
                        value = 0.5 * _vec1 ( icmp * numTotalDof1 + iDofEdge ) + _vec2[iDof2];
                        //                                            std::cout << " transfering " << value << " from P2 " << iDofEdge << " to P1 " << iDof2 << " ( " << iter->second << " ) " << std::endl;
                        _vec2.setCoefficient ( iDof2, value );
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
        for ( ID iEdge = 0; iEdge < nBEdges2; ++iEdge )
        {
            if ( (int) _fespace2.mesh()->edgeList (iEdge).markerID() != M_data->fluidInterfaceFlag() )
            {
                continue;
            }
            // Now that we have an edge on the FS interface, let's get its ID
            ID edgeID = _fespace2.mesh()->edgeList (iEdge).id();
            // dof position of the edge since unknowns are store using [ node | edges | faces | volumes ]
            int iDofEdge = numTotalVert2 + edgeID;

            for ( ID l = 0; l < nDofPerEdge1; ++l )
            {
                // Loop on data vector components
                for ( UInt icmp = 0; icmp < nDimensions; ++icmp )
                {
                    // interpolation of the nodal values in the middle of the segement
                    ID node1 = _fespace2.mesh()->edgeList (iEdge).point (0).id();
                    ID node2 = _fespace2.mesh()->edgeList (iEdge).point (1).id();
                    value = 0.5 * (_vec2 (icmp * numTotalDof2 + node1) + _vec2 (icmp * numTotalDof2 + node2) );
                    _vec2.setCoefficient (iDofEdge, value);
                }
            }
        }
    }


}


} // Namespace LifeV

