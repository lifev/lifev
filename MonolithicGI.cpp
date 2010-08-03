/* -*- mode: c++ -*-

   This file is part of the LifeV library

   Author(s): Paolo Crosetto <crosetto@iacspc70.epfl.ch>
   Date: 2008-09-17

   Copyright (C) 2008

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file (monolithicGI)
   \author crosetto <Paolo Crosetto>
   \date (17/09/2008)
*/

#include <lifemc/lifesolver/MonolithicGI.hpp>

namespace LifeV
{

MonolithicGI::MonolithicGI():
    super(),
    M_mapWithoutMesh(),
    M_domainVelImplicit(true),
    M_convectiveTermDer(true),
    M_solidOper(),
    M_fluidOper(),
    M_meshOper(),
    M_interface(0),
    M_compPrecPtr(),
    M_precMatrPtr()
{}

void
MonolithicGI::setupFluidSolid()
{
    super::setupFluidSolid();
    UInt offset = M_monolithicMap->getMap(Unique)->NumGlobalElements();
    M_mapWithoutMesh.reset(new EpetraMap(*M_monolithicMap));
    M_BCh_mesh->setOffset(offset);

    *this->M_monolithicMap += this->M_mmFESpace->map();
    /* OBSOLETE
       if(M_data->dataFluid()->useShapeDerivatives())
       {
       M_epetraOper.reset( new Epetra_FullMonolithic(this));
       M_solid->setOperator(*M_epetraOper);
       }*/
    //std::cout<<"map global elements : "<<M_monolithicMap->getMap(Unique)->NumGlobalElements()<<std::endl;
    M_interface=M_monolithicMatrix->getInterface();
    M_monolithicInterfaceMap=M_monolithicMatrix->getInterfaceMap();

    vector_type u0(*this->M_monolithicMap);
    M_bdf.reset(new BdfT<vector_type>(M_data->dataFluid()->dataTime()->getBDF_order()));
    M_bdf->initialize_unk(u0);
    this->M_rhs.reset(new vector_type(*this->M_monolithicMap));
    this->M_rhsFull.reset(new vector_type(*this->M_monolithicMap));
    M_uk.reset (new vector_type(*this->M_monolithicMap));
    M_un.reset (new vector_type(*this->M_monolithicMap));

    M_meshMotion.reset(new meshmotion_raw_type(*M_mmFESpace,
                                               *M_epetraComm,
                                               *M_monolithicMap,
                                               offset));
    M_fluid.reset     (new fluid_raw_type(dataFluid(),
                                          *M_uFESpace,
                                          *M_pFESpace,
                                          *M_mmFESpace,
                                          *M_epetraComm,
                                          *M_monolithicMap));
    M_solid.reset     (new solid_raw_type(dataSolid(),
                                          *M_dFESpace,
                                          *M_epetraComm,
                                          *M_monolithicMap,
                                          M_offset
                                          ));
}

void
MonolithicGI::updateSystem()
{
    //M_meshMotion->dispOld() is at time n-1 !!
    UInt offset(M_solidAndFluidDim + nDimensions*M_interface);
    vector_ptrtype meshDispDiff(new vector_type(M_mmFESpace->map()));
    meshDispDiff->subset(*M_uk, offset); //if the conv. term is to be condidered implicitly
    M_meshMotion->initialize(*meshDispDiff);//M_disp is set to the total mesh disp.`
    super::updateSystem();
    M_un.reset(new vector_type(*M_uk));
}

void
MonolithicGI::couplingMatrix(matrix_ptrtype &fullMMatrix, const std::map<ID, ID>& locDofMap, int coupling)
{
    if(coupling-16>=0)
    {
        UInt solidFluidInterface(M_solidAndFluidDim + nDimensions*M_interface);
        std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
        std::map<ID, ID>::const_iterator ITrow;
        for(UInt dim = 0; dim < nDimensions; ++dim)
        {
            for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
            {
                if(M_interfaceMap->getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
                    fullMMatrix->set_mat_inc(solidFluidInterface + ITrow->first + dim*M_mmFESpace->dof().numTotalDof() - 1, M_offset + ITrow->second-1 + dim* M_dFESpace->dof().numTotalDof(), (-1.0)*M_data->dataFluid()->dataTime()->getTimeStep()*M_solid->rescaleFactor()/**1.e-2*//*scaling of the solid matrix*/ );
                }
            }
        }
        coupling-=16;
    }
    couplingMatrix2(fullMMatrix, M_dofStructureToHarmonicExtension->locDofMap(), coupling);
}

void
MonolithicGI::couplingMatrix2(matrix_ptrtype &bigMatrix, const std::map<ID, ID>& locDofMap, int coupling)
{
    std::map<ID, ID>::const_iterator ITrow;
    int flag(coupling);

    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow, flag=coupling)
        {
            if(M_interfaceMap->getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                if(flag-8>=0)//right low
                {
                    bigMatrix->set_mat_inc( M_offset + ITrow->second-1 + dim* M_dFESpace->dof().numTotalDof(),(int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (1.) );//right low
                    flag -= 8;
                }
                if(flag-4>=0)// right up
                {
                    bigMatrix->set_mat_inc( ITrow->first-1 + dim* M_uFESpace->dof().numTotalDof(), (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (-1.) );//right up
                    flag -= 4;
                }
                if(flag-2>=0)//low left
                {
                    bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), 1.);//low left
                    flag -= 2;
                }
                if(flag-1>=0)//low right
                    bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (M_offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), (-1.*M_solid->rescaleFactor()/*/M_data->dataFluid()->timestep()*/));//low right

                bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim , (int)(*this->M_numerationInterface)[ITrow->second /*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, 0.0);
            }
        }
    }


}

void
MonolithicGI::buildSystem ()
{
    super::buildSystem();
    M_meshMotion->computeMatrix();
}

void
MonolithicGI::evalResidual( vector_type&       res,
                              const vector_type& disp,
                              const UInt          iter )
{
    if(iter == 0)
    {
        M_nbEval = 0; // new time step
    }
    else
    {
        this->M_solid->updateVel();
    }
    M_uk.reset(new vector_type(disp));
    this->M_beta.reset(new vector_type(M_uFESpace->map()));
    UInt offset(M_solidAndFluidDim + nDimensions*M_interface);

    vector_ptrtype meshDispDiff(new vector_type(M_mmFESpace->map()));
    vector_ptrtype meshDispOld(new vector_type(M_mmFESpace->map()));

    meshDispDiff->subset(disp, offset); //if the conv. term is to be condidered implicitly
    meshDispOld->subset(*M_un, offset);
    //meshDispDiff->subset(*M_uk, offset); //if the mesh motion is at the previous nonlinear step (FP) in the convective term
    //meshDispDiff->subset(*M_un, offset); //if we linearize in a semi-implicit way
    M_meshMotion->initialize(*meshDispDiff);//M_disp is set to the total mesh disp.
    double alpha = 1/M_data->dataFluid()->dataTime()->getTimeStep();
    vector_type mmRep(*meshDispDiff, Repeated);// just to repeat dispDiff. No way witout copying?
    this->moveMesh(mmRep);// re-initialize the mesh points
    *meshDispDiff -= *meshDispOld;//relative displacement
    if(!M_domainVelImplicit)
    {
        meshDispDiff=meshDispOld;// at time n /*->subset(*M_un, offset)*/; //if the mesh motion is at the previous time step in the convective term
        *meshDispDiff -= M_meshMotion->dispOld();//at time n-1
    }
    *meshDispDiff *= -alpha;// -w, mesh velocity
    mmRep = *meshDispDiff;

    this->interpolateVelocity(mmRep, *this->M_beta);
    //            *this->M_beta *= -alpha; //HE solution scaled!
    vector_ptrtype fluid(new vector_type(this->M_uFESpace->map()));
    if(!M_convectiveTermDer)
        fluid->subset(*M_un/**M_unOld*/, 0);
    else
        fluid->subset(disp, 0);
    *this->M_beta += *fluid/*M_un or disp, it could be also M_uk in a FP strategy*/;

    assembleFluidBlock(iter);
    assembleMeshBlock(iter);

        if(!M_monolithicMatrix->set())
        {
            M_BChs.push_back(M_BCh_d);
            M_BChs.push_back(M_BCh_u);
            M_FESpaces.push_back(M_dFESpace);
            M_FESpaces.push_back(M_uFESpace);

            M_BChs.push_back(M_BCh_mesh);
            M_FESpaces.push_back(M_mmFESpace);

            M_monolithicMatrix->push_back_matrix(M_solidBlock, false);
            M_monolithicMatrix->push_back_matrix(M_fluidBlock, true);
            M_monolithicMatrix->push_back_matrix(M_meshBlock, false);
            M_monolithicMatrix->setConditions(M_BChs);
            M_monolithicMatrix->setSpaces(M_FESpaces);
            M_monolithicMatrix->setOffsets(3, M_offset, 0, M_solidAndFluidDim + nDimensions*M_interface);
            M_monolithicMatrix->coupler(M_monolithicMap, M_dofStructureToHarmonicExtension->locDofMap(), M_numerationInterface, M_data->dataFluid()->dataTime()->getTimeStep());
        }
        else
        {
            M_monolithicMatrix->replace_matrix(M_solidBlock, 0);
            M_monolithicMatrix->replace_matrix(M_fluidBlock, 1);
            M_monolithicMatrix->replace_matrix(M_meshBlock, 2);
        }

        M_monolithicMatrix->blockAssembling();

        if ( !M_BCh_u->bdUpdateDone() )
            M_BCh_u->bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );

        M_BCh_d->setOffset(M_offset);
        if ( !M_BCh_d->bdUpdateDone() )
            M_BCh_d->bdUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
        M_BCh_mesh->setOffset(M_solidAndFluidDim + nDimensions*M_interface);
        if ( !M_BCh_mesh->bdUpdateDone() )
            M_BCh_mesh->bdUpdate( *M_mmFESpace->mesh(), M_mmFESpace->feBd(), M_mmFESpace->dof() );


        M_monolithicMatrix->applyBoundaryConditions(dataFluid().dataTime()->getTime(), M_rhsFull);
        M_monolithicMatrix->GlobalAssemble();


    if(M_data->DDBlockPreconditioner()>=2 && M_data->DDBlockPreconditioner()!=5 && M_data->DDBlockPreconditioner()!=7 && M_data->DDBlockPreconditioner()!=8 && M_data->DDBlockPreconditioner()<11)
    {
        if(!M_robinCoupling.get())
        {
            M_robinCoupling.reset( new matrix_type(*M_monolithicMap));//this matrix is a preconditioner that mixes the coupling conditions at the interface
            M_robinCoupling->insertValueDiagonal(1., *M_monolithicMap);
            robinCoupling(M_robinCoupling, M_data->RobinNeumannFluidCoefficient(), M_data->RobinNeumannSolidCoefficient());
            M_robinCoupling->GlobalAssemble();
        }
    }


    if(M_data->DDBlockPreconditioner()!=5 && M_data->DDBlockPreconditioner()<9)
        M_precMatrPtr.reset(new matrix_type(*M_monolithicMap));
    super::evalResidual( disp, M_rhsFull, res, false );
}

int MonolithicGI::setupBlockPrec(vector_type& rhs)
{
    boost::shared_ptr<IfpackComposedPrec>  ifpackCompPrec;

    if(M_data->dataFluid()->useShapeDerivatives())
    {
        if(M_data->DDBlockPreconditioner()==6)
        {
            M_precMatrPtr.reset(new matrix_type(*M_monolithicMatrix->getMatrix()));
        }
        matrix_ptrtype monolithicMatrix(new matrix_type(*M_monolithicMap));
        *monolithicMatrix+=*M_monolithicMatrix->getMatrix();
        shapeDerivatives(monolithicMatrix,*M_uk /*subX*/, M_domainVelImplicit, M_convectiveTermDer);
        monolithicMatrix->GlobalAssemble();
        M_monolithicMatrix->getMatrix()=monolithicMatrix;
    }
    switch(M_data->DDBlockPreconditioner())
    {
    case 1:
        {
            M_meshMotion->setMatrix(M_precMatrPtr);
        }
        break;
    case 2:    case 3:    case 4:    case 5:
        {
        }
        break;
    case 6:
        {
            M_precMatrPtr->GlobalAssemble();
        }
        break;
    case 7:
    case 8:
        {
            displayer().leaderPrint("Preconditioner type not yet implemented," );
            displayer().leaderPrint("change the entry DDBlockPrec in the data file");
            throw WRONG_PREC_EXCEPTION();
        }
        break;
    case 9://like 6 but with composed prec.
        {
            M_fluidBlock.reset(new matrix_type(*M_monolithicMap));
            M_fluid->getFluidMatrix( *M_fluidBlock);
            M_fluid->updateStab( *M_fluidBlock);//applies the stabilization terms
            *M_fluidBlock+=*M_solidBlock;
            M_fluidBlock->insertValueDiagonal(1.,  M_mmFESpace->map(), mapWithoutMesh().getMap(Unique)->NumGlobalElements());
            couplingMatrix(M_fluidBlock, M_dofStructureToHarmonicExtension->locDofMap(), 15);

            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_u, M_uFESpace->feBd(), 1., dataSolid().dataTime()->getTime() );
            M_fluidBlock->GlobalAssemble();

            M_BCh_d->setOffset(M_offset);
            bcManageMatrix( *M_fluidBlock, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataSolid().dataTime()->getTime() );
            //super::applyPreconditioner(M_robinCoupling, rhs);
            applyPreconditioner(M_robinCoupling, M_fluidBlock);
            M_fluidBlock->GlobalAssemble();
            M_fluidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_fluidBlock));

            ifpackCompPrec =  boost::dynamic_pointer_cast< IfpackComposedPrec, prec_raw_type > (M_compPrecPtr);

            if(ifpackCompPrec->set() && ifpackCompPrec->getNumber())
            {
                ifpackCompPrec->replace(M_fluidOper, 0);
            }
            else
            {
                M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap));
                *M_solidBlockPrec += *M_meshBlock;
                M_solidBlockPrec->insertValueDiagonal(1., mapWithoutMesh(), 0);
                couplingMatrix(M_solidBlockPrec, M_dofStructureToHarmonicExtension->locDofMap(), 16);
                M_solidBlockPrec->GlobalAssemble();
                M_solidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_solidBlockPrec));
                ifpackCompPrec->buildPreconditioner(M_fluidOper);
                ifpackCompPrec->push_back(M_solidOper);
            }
        }
        break;
    case 10:
        {
            M_fluidBlock.reset(new matrix_type(*M_monolithicMap));
            M_fluid->getFluidMatrix( *M_fluidBlock);
            M_fluid->updateStab( *M_fluidBlock);//applies the stabilization terms
            *M_fluidBlock+=*M_solidBlock;
             M_fluidBlock->insertValueDiagonal(1., M_mmFESpace->map(), mapWithoutMesh().getMap(Unique)->NumGlobalElements());
            couplingMatrix(M_fluidBlock, M_dofStructureToHarmonicExtension->locDofMap(), 15);
            //*M_fluidBlock+=*M_matrShapeDer;
            shapeDerivatives(M_fluidBlock,*M_uk /*subX*/, M_domainVelImplicit, M_convectiveTermDer);

            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_u, M_uFESpace->feBd(), 1., dataSolid().dataTime()->getTime() );

            M_BCh_d->setOffset(M_offset);
            bcManageMatrix( *M_fluidBlock, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataSolid().dataTime()->getTime() );

            //super::applyPreconditioner(M_robinCoupling, rhs);
            applyPreconditioner(M_robinCoupling, M_fluidBlock);
            M_fluidBlock->GlobalAssemble();
            M_fluidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_fluidBlock));

            ifpackCompPrec =  boost::dynamic_pointer_cast< IfpackComposedPrec, prec_raw_type > (M_compPrecPtr);

            if(ifpackCompPrec->set() && ifpackCompPrec->getNumber())
            {
                ifpackCompPrec->replace(M_fluidOper, 1);
            }
            else
            {
                M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap));
                *M_solidBlockPrec += *M_meshBlock;
                M_solidBlockPrec->insertValueDiagonal(1., mapWithoutMesh(), 0);
                //M_meshMotion->setMatrix(M_solidBlockPrec);
                M_solidBlockPrec->GlobalAssemble();

                M_solidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_solidBlockPrec));

                ifpackCompPrec->buildPreconditioner(M_solidOper);
                ifpackCompPrec->push_back(M_fluidOper);
            }
        }
	    break;
    case 11://P1, P2, P3
        {
            M_fluidBlock.reset(new matrix_type(*M_monolithicMap));
            M_fluid->getFluidMatrix( *M_fluidBlock);
            M_fluid->updateStab( *M_fluidBlock);//applies the stabilization terms
            M_fluidBlock->insertValueDiagonal(1., M_mmFESpace->map(), mapWithoutMesh().getMap(Unique)->NumGlobalElements());
            M_fluidBlock->insertValueDiagonal(1., M_dFESpace->map(), M_offset);
            couplingMatrix(M_fluidBlock, M_dofStructureToHarmonicExtension->locDofMap(), 7);
            //*M_fluidBlock+=*M_matrShapeDer;
            shapeDerivatives(M_fluidBlock,*M_uk /*subX*/, M_domainVelImplicit, M_convectiveTermDer);

            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_u, M_uFESpace->feBd(), 1., dataSolid().dataTime()->getTime() );

            M_fluidBlock->GlobalAssemble();

            M_fluidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_fluidBlock));

            ifpackCompPrec =  boost::dynamic_pointer_cast< IfpackComposedPrec, prec_raw_type > (M_compPrecPtr);

            if(ifpackCompPrec->set() && ifpackCompPrec->getNumber())
            {
                ifpackCompPrec->replace(M_fluidOper, 2);
            }
            else
            {
                matrix_ptrtype mesh(new matrix_type(*M_monolithicMap));
                *mesh += *M_meshBlock;
                mesh->insertValueDiagonal(1., mapWithoutMesh(), 0);
                //M_meshMotion->setMatrix(M_meshBlock);
                couplingMatrix(mesh, M_dofStructureToHarmonicExtension->locDofMap(), 16);
                mesh->GlobalAssemble();


                M_meshOper.reset(new IfpackComposedPrec::operator_raw_type(*mesh));


                M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap, 1));
                *M_solidBlockPrec += *M_solidBlock;
                M_solidBlockPrec->insertValueDiagonal(1., M_uFESpace->map() );
                M_solidBlockPrec->insertValueDiagonal(1., M_pFESpace->map()+M_fluxes , M_uFESpace->dof().numTotalDof()*nDimensions );
                M_solidBlockPrec->insertValueDiagonal(1., *M_monolithicInterfaceMap, M_solidAndFluidDim);
                M_solidBlockPrec->insertValueDiagonal(1.,  M_mmFESpace->map(), M_solidAndFluidDim + nDimensions*M_interface);
//                 bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataSolid().dataTime()->getTime() );
                M_BCh_d->setOffset(M_offset);
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataSolid().dataTime()->getTime() );
                M_solidBlockPrec->GlobalAssemble();

                M_solidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_solidBlockPrec));


                ifpackCompPrec->buildPreconditioner(M_solidOper);
                ifpackCompPrec->push_back(M_meshOper);
                ifpackCompPrec->push_back(M_fluidOper);
            }
        }
        break;
    case 12://P1, P2, P3
        {
            M_fluidBlock.reset(new matrix_type(*M_monolithicMap));
            M_fluid->getFluidMatrix( *M_fluidBlock);
            M_fluid->updateStab( *M_fluidBlock);//applies the stabilization terms
            M_fluidBlock->insertValueDiagonal(1., M_mmFESpace->map(), mapWithoutMesh().getMap(Unique)->NumGlobalElements());
            M_fluidBlock->insertValueDiagonal(1., M_dFESpace->map(), M_offset);
            couplingMatrix(M_fluidBlock, M_dofStructureToHarmonicExtension->locDofMap(), 7);
            //this->shapeDerivatives(M_fluidBlock,*M_uk /*subX*/);

            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_u, M_uFESpace->feBd(), 1., dataSolid().dataTime()->getTime() );
            M_fluidBlock->GlobalAssemble();
            M_fluidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_fluidBlock));

            ifpackCompPrec =  boost::dynamic_pointer_cast< IfpackComposedPrec, prec_raw_type > (M_compPrecPtr);

            if( ifpackCompPrec->set() && ifpackCompPrec->getNumber())
            {
                ifpackCompPrec->replace(M_fluidOper, 2);
            }
            else
            {
                M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap, 1));
                *M_solidBlockPrec += *M_solidBlock;
                M_solidBlockPrec->insertValueDiagonal(1., M_uFESpace->map() );
                M_solidBlockPrec->insertValueDiagonal(1., M_pFESpace->map(), M_uFESpace->dof().numTotalDof()*nDimensions );
                M_solidBlockPrec->insertValueDiagonal(1., *M_monolithicInterfaceMap, M_solidAndFluidDim);
                M_solidBlockPrec->insertValueDiagonal(1., M_mmFESpace->map(), M_solidAndFluidDim + nDimensions*M_interface);
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataSolid().dataTime()->getTime() );
                M_solidBlockPrec->GlobalAssemble();

                M_solidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_solidBlockPrec));

                matrix_ptrtype mesh(new matrix_type(*M_monolithicMap));
                *mesh += *M_meshBlock;
                mesh->insertValueDiagonal(1., mapWithoutMesh(), 0);
                couplingMatrix(mesh, M_dofStructureToHarmonicExtension->locDofMap(), 16);
                mesh->GlobalAssemble();
                shapeDerivatives(mesh,*M_uk /*subX*/, M_domainVelImplicit, M_convectiveTermDer);
                mesh->GlobalAssemble();

                M_meshOper.reset(new IfpackComposedPrec::operator_raw_type(*mesh));

                ifpackCompPrec->buildPreconditioner(M_meshOper);
                ifpackCompPrec->push_back(M_solidOper);
                ifpackCompPrec->push_back(M_fluidOper);
            }
        }
        break;
    default:
        {}
        break;
    }
}

void MonolithicGI::solveJac(vector_type       &_step,
                              const vector_type &_res,
                              const Real       _linearRelTol)
{
    setMatrix(*M_monolithicMatrix->getMatrix());
    setupBlockPrec(*const_cast<vector_type*>(&_res));


    M_solid->getDisplayer().leaderPrint("  M-  Jacobian NormInf res:                    ", _res.NormInf(), "\n");
    M_solid->getDisplayer().leaderPrint("  M-  Solving Jacobian system ...              \n" );

    switch(M_data->DDBlockPreconditioner())
    {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
        M_linearSolver->buildPreconditioner(M_precMatrPtr);
        this->iterateFullMonolithic(_res, _step, M_precMatrPtr->getMatrixPtr(), M_linearSolver);
        break;

    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
        this->iterateFullMonolithic(_res, _step, M_compPrecPtr->getPrecPtr(),     M_linearSolver);
        break;

    default:
        if(!M_precMatrPtr.get())
        {
            displayer().leaderPrint("Preconditioner type not implemented," );
            displayer().leaderPrint("change the entry DDBlockPrec in the data file");
            throw WRONG_PREC_EXCEPTION();
        }
        else
        {
            displayer().leaderPrint("default (Amesos) preconditioner choice" );
        }
        this->iterateFullMonolithic(_res, _step, M_precMatrPtr->getMatrixPtr(), M_linearSolver);
        break;
    }

    M_solid->getDisplayer().leaderPrint("  M-  Jacobian NormInf res:                    ", _step.NormInf(), "\n");

}



void MonolithicGI::initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                                 FSIOperator::solid_type::value_type::Function const& p0,
                                 FSIOperator::solid_type::value_type::Function const& d0,
                                 FSIOperator::solid_type::value_type::Function const& w0,
                                 FSIOperator::solid_type::value_type::Function const& df0 )
{
    super::initialize(u0, p0, d0, w0, df0);

    vector_type df(M_mmFESpace->map());
    M_mmFESpace->interpolate(df0, df, M_data->dataSolid()->dataTime()->getTime());

    M_un->add(df, M_solidAndFluidDim+getDimInterface());

    M_meshMotion->setDisplacement(df);
}

void MonolithicGI::shapeDerivatives(matrix_ptrtype sdMatrix, const vector_type& sol, bool domainVelImplicit, bool convectiveTermDer)
{
    double alpha = 1./M_data->dataFluid()->dataTime()->getTimeStep();
    vector_ptrtype rhsNew(new vector_type(*M_monolithicMap));
    vector_type un(M_uFESpace->map()/*+M_pFESpace->map()*/);
    vector_type uk(M_uFESpace->map()+M_pFESpace->map());

    //vector_type meshVel(M_meshMotion->dispDiff(), Repeated);
    vector_ptrtype meshVel;//(M_mmFESpace->map(), Repeated);
    meshVel.reset(new vector_type(M_mmFESpace->map()));

    UInt offset(M_solidAndFluidDim + nDimensions*M_interface);
    if(domainVelImplicit)
    {
        vector_type meshDispOld(M_mmFESpace->map());
        meshVel->subset(sol, offset); //if the conv. term is to be condidered implicitly
        meshDispOld.subset(*M_un, offset);
        *meshVel -= meshDispOld;
    }
    else
    {
        meshVel->subset(*M_un, offset); //if the conv. term is to be condidered partly explicitly
        *meshVel -= M_meshMotion->dispOld();
    }

    if(convectiveTermDer)
        un.subset(sol, 0);
    else
        un.subset(*this->M_un, 0);


    *meshVel *= alpha;
    vector_ptrtype meshVelRep(new vector_type(M_mmFESpace->map(), Repeated));
    *meshVelRep = *meshVel;

    uk.subset(sol, 0);
    vector_type dvfm(M_uFESpace->map(), Repeated);
    vector_type vfm(M_uFESpace->map(), Repeated);
    this->transferMeshMotionOnFluid(*meshVelRep, vfm);

    M_fluid->updateShapeDerivatives(*sdMatrix,
                                    alpha,
                                    un,//un if !domainVelImplicit, otherwise uk
                                    uk,//uk
                                    vfm, //(xk-xn)/dt (FI), or (xn-xn-1)/dt (CE)//Repeated
                                    M_solidAndFluidDim+M_interface*nDimensions,
                                    *M_uFESpace,
                                    domainVelImplicit,
                                    convectiveTermDer
                                    );
}

void
MonolithicGI::setUp( const GetPot& dataFile )
{
    super::setUp(dataFile);

    M_domainVelImplicit     = dataFile( "fluid/domainVelImplicit", true);
    M_convectiveTermDer     = dataFile( "fluid/convectiveTermDer", false);

    M_compPrecPtr.reset(new IfpackComposedPrec(M_epetraComm.get()));
    M_compPrecPtr->setDataFromGetPot(dataFile, "linear_system/prec");//to avoid if we build the prec from a matrix.
}


void MonolithicGI::applyPreconditioner( const matrix_ptrtype prec, matrix_ptrtype& oper )
{
    matrix_type tmpMatrix(prec->getMap(), 1);
    EpetraExt::MatrixMatrix::Multiply( *prec->getMatrixPtr(),
                                       false,
                                       *oper->getMatrixPtr(),
                                       false,
                                       *tmpMatrix.getMatrixPtr());
    oper->swapCrsMatrix(tmpMatrix);
}

void
MonolithicGI::assembleFluidBlock(UInt iter)
{
    double alpha = 1./M_data->dataFluid()->dataTime()->getTimeStep();
    M_fluid->recomputeMatrix(true);
    M_fluid->updateSystem(alpha, *this->M_beta, *this->M_rhs );//here it assembles the fluid matrices
    if(iter==0)
    {
        *this->M_rhs += this->M_fluid->matrMass()*this->M_bdf->time_der( M_data->dataFluid()->dataTime()->getTimeStep() );// fluid time discr.
        super::couplingRhs( this->M_rhs, this->M_un);//adds to the solid rhs the terms due to the coupling
        super::updateSolidSystem(this->M_rhs);//updates the rhs in the solid system (the solid time discr. is implemented there)
    }
    M_fluidBlock.reset(new matrix_type(*M_monolithicMap));
    M_fluid->getFluidMatrix( *M_fluidBlock);
    M_fluid->updateStab( *M_fluidBlock);//applies the stabilization terms
    *this->M_rhsFull = *this->M_rhs;//this will be the rhs with the BCs (the one who changes at each nonlinear iteration)
}

void
MonolithicGI::assembleMeshBlock(UInt iter)
{
    M_meshBlock.reset(new matrix_type(*M_monolithicMap));
    M_meshMotion->setMatrix(M_meshBlock);
    M_meshBlock->GlobalAssemble();
    UInt offset(M_solidAndFluidDim+nDimensions*M_interface);
    std::map<ID, ID>::const_iterator ITrow;
    std::map<ID, ID> locdofmap(M_dofStructureToHarmonicExtension->locDofMap());

    for ( ID dim=0; dim < nDimensions; ++dim )
        for( ITrow = locdofmap.begin(); ITrow != locdofmap.end(); ++ITrow )
        {
            UInt i = ITrow->first;
            M_meshBlock->diagonalize(i+offset+dim*M_mmFESpace->dof().numTotalDof()-1 , 1.);
        }
}


 void
MonolithicGI::robinCoupling( matrix_ptrtype& IdentityMatrix, const Real& alphaf, const Real& alphas, int coupling ) // not working with non-matching grids
{
    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    //Real* entry(new Real(1.));

    //addDiagonalEntries(alphaf, IdentityMatrix, *M_solidInterfaceMap, M_offset, true);
    if(coupling-4 >= 0)
    {
            for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
            {
                if(M_interfaceMap->getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
                    for(UInt dim = 0; dim < nDimensions; ++dim)
                    {
                        IdentityMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second ] - 1 + dim*M_interface + M_solidAndFluidDim, (M_offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), alphas);//low right//to enable
                    }
                }
            }
            coupling -= 4;
    }
    if(coupling-2 >= 0)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap->getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
            {
                for(UInt dim = 0; dim < nDimensions; ++dim)
                {
                    IdentityMatrix->set_mat_inc( ITrow->first-1 + dim* M_uFESpace->dof().numTotalDof(), (int)(*this->M_numerationInterface)[ITrow->second ] - 1 + dim*M_interface + M_solidAndFluidDim, alphaf );//right up
                }
            }
        }
        coupling -= 2;
    }
    if(coupling-1 >= 0)
    {
            for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
            {
                if(M_interfaceMap->getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
                    for(UInt dim = 0; dim < nDimensions; ++dim)
                    {
                    IdentityMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second ] - 1 + dim*M_interface + M_solidAndFluidDim, (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), alphas);//low left
                    }
                }
            }
    }
}

namespace
{
BlockMatrix*    createAdditiveSchwarzGI(){ return new BlockMatrix(31); }
FSIOperator*    createFM(){ return new MonolithicGI(); }
}

bool MonolithicGI::reg =  BlockPrecFactory::instance().registerProduct("AdditiveSchwarzGI"  , &createAdditiveSchwarzGI ) &&
    BlockMatrix::Factory::instance().registerProduct("AdditiveSchwarzGI"  , &createAdditiveSchwarzGI )
&&
    FSIFactory::instance().registerProduct( "monolithicGI", &createFM );
}
