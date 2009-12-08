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
   \file (fullMonolithic)
   \author crosetto <Paolo Crosetto>
   \date (17/09/2008)
*/

#include <lifemc/lifesolver/fullMonolithic.hpp>

namespace LifeV
{

fullMonolithic::fullMonolithic():
    super(),
    M_mapWithoutMesh(),
    M_domainVelImplicit(true),
    M_convectiveTermDer(true)
{}

void
fullMonolithic::setUp( const GetPot& dataFile )
{
    super::setUp(dataFile);
    M_domainVelImplicit     = dataFile( "fluid/domainVelImplicit", true);
    M_convectiveTermDer     = dataFile( "fluid/convectiveTermDer", false);
}

void
fullMonolithic::setupFEspace()
{
	super::setupFEspace();
}

void
fullMonolithic::setupFluidSolid()
{
    super::setupFluidSolid();
    UInt offset = M_monolithicMap->getMap(Unique)->NumGlobalElements();
    M_mapWithoutMesh.reset(new EpetraMap(*M_monolithicMap));
    M_BCh_mesh->setOffset(offset);

    *this->M_monolithicMap += this->M_mmFESpace->map();
    /* OBSOLETE
       if(M_dataFluid->useShapeDerivatives())
       {
       M_epetraOper.reset( new Epetra_FullMonolithic(this));
       M_solid->setOperator(*M_epetraOper);
       }*/
    //std::cout<<"map global elements : "<<M_monolithicMap->getMap(Unique)->NumGlobalElements()<<std::endl;
    vector_type u0(*this->M_monolithicMap);
    M_bdf.reset(new BdfT<vector_type>(M_dataFluid->getBDF_order()));
    M_bdf->initialize_unk(u0);
    this->M_rhs.reset(new vector_type(*this->M_monolithicMap));
    this->M_rhsFull.reset(new vector_type(*this->M_monolithicMap));
    M_un.reset (new vector_type(*this->M_monolithicMap));
    M_meshMotion.reset(new FSIOperator::meshmotion_raw_type(*M_mmFESpace,
                                                            *M_epetraComm,
                                                            *M_monolithicMap,
                                                            offset));
    M_fluid.reset(new FSIOperator::fluid_raw_type(dataFluid(),
                                                  *M_uFESpace,
                                                  *M_pFESpace,
                                                  *M_mmFESpace,
                                                  *M_epetraComm,
                                                  *M_monolithicMap));
    M_solid.reset(new solid_raw_type(dataSolid(),
                                     *M_dFESpace,
                                     *M_epetraComm,
                                     *M_monolithicMap,
                                     M_offset
                                     ));

}

void
fullMonolithic::updateSystem(const vector_type& solution)
{
    //M_meshMotion->dispOld() is at time n-1 !!
    UInt offset(M_solidAndFluidDim + nDimensions*M_interface);
    vector_ptrtype meshDispDiff(new vector_type(M_mmFESpace->map()));
    meshDispDiff->subset(*M_un, offset); //if the conv. term is to be condidered implicitly
    M_meshMotion->initialize(*meshDispDiff);//M_disp is set to the total mesh disp.`
    super::updateSystem(solution);
}

void
fullMonolithic::couplingMatrix(matrix_ptrtype         &fullMMatrix, int coupling)
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
                if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
                    fullMMatrix->set_mat_inc(solidFluidInterface + ITrow->first + dim*M_mmFESpace->dof().numTotalDof() - 1, M_offset + ITrow->second-1 + dim* M_dFESpace->dof().numTotalDof(), (-1.0)*M_dataFluid->getTimeStep()*M_solid->rescaleFactor()/**1.e-2*//*scaling of the solid matrix*/ );
                }
            }
        }
        coupling-=16;
    }
    super::couplingMatrix(fullMMatrix, coupling);
}

void
fullMonolithic::buildSystem ()
{
    super::buildSystem();
    M_meshMotion->computeMatrix();
    //    double dt = M_dataFluid->timestep();
    //    M_meshMotion->rescaleMatrix(dt);
}
void
fullMonolithic::evalResidual( vector_type&       res,
                              const vector_type& disp,
                              const UInt          iter )
{
    if(iter == 0)
    {
        M_nbEval = 0; // new time step
    }
    else
    {
        this->M_solid->updateStuff();
    }
    if((iter==0)|| !this->M_dataFluid->isSemiImplicit())
    {
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
        double alpha = 1/M_dataFluid->getTimeStep();
        vector_type mmRep(*meshDispDiff, Repeated);// just to repeat dispDiff. No way witout copying?
        this->moveMesh(mmRep);// re-initialize the mesh points
        *meshDispDiff -= *meshDispOld;//relative displacement
        if(!M_domainVelImplicit)
        {
            meshDispDiff=meshDispOld;// at time n /*->subset(*M_un, offset)*/; //if the mesh motion is at the previous time step in the convective term
            //        std::cout<<"meshDispDiff1.5 "<<meshDispDiff->NormInf()<<std::endl;
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
        //          if(firstIter)
        M_fluid->recomputeMatrix(true);
        M_fluid->updateSystem(alpha, *this->M_beta, *this->M_rhs );//here it assembles the fluid matrices
        if(iter==0)
        {
            *this->M_rhs += this->M_fluid->matrMass()*this->M_bdf->time_der( M_dataFluid->getTimeStep() );// fluid time discr.
            super::couplingRhs( this->M_rhs, this->M_un);//adds to the solid rhs the terms due to the coupling
            super::updateSolidSystem(this->M_rhs);//updates the rhs in the solid system (the solid time discr. is implemented there)
        }
        M_monolithicMatrix.reset(new matrix_type(*M_monolithicMap/*, this->M_fluid->getMeanNumEntries()*/));
        M_fluid->getFluidMatrix( *M_monolithicMatrix);// adds the fluid part
        M_fluid->updateStab( *M_monolithicMatrix);//adds the stabilization terms
        this->updateMatrix(*M_monolithicMatrix);// adds the solid and the coupling parts
        this->M_rhsFull.reset(new vector_type(*this->M_rhs));//this will be the rhs with the BCs (the one who changes at each nonlinear iteration)
        M_meshMotion->applyBoundaryConditions(*M_rhsFull, *M_BCh_mesh);
        M_meshMotion->setMatrix(M_monolithicMatrix);// adds the mesh part


	    M_BCh_flux->setOffset(M_offset-M_fluxes);
	    if ( !M_BCh_flux->bdUpdateDone() )
            M_BCh_flux->bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );
	    bcManage( *M_monolithicMatrix, *this->M_rhsFull, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_flux, M_uFESpace->feBd(), 1., dataSolid().getTime() );

	    M_BCh_Robin->setOffset(M_offset);
	    if ( !M_BCh_Robin->bdUpdateDone() )
            M_BCh_Robin->bdUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );
	    bcManage( *M_monolithicMatrix, *this->M_rhsFull, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataSolid().getTime() );

        super::evalResidual( *M_BCh_u, *M_BCh_d, disp, this->M_rhsFull, res, false);
        if(M_DDBlockPrec>=2 && M_DDBlockPrec!=5 && M_DDBlockPrec!=7 && M_DDBlockPrec!=8 && M_DDBlockPrec<11)
        {
            if(!M_robinCoupling.get())
            {
                M_robinCoupling.reset( new matrix_type(*M_monolithicMap));//this matrix is a preconditioner that mixes the coupling conditions at the interface
                super::robinCoupling(M_robinCoupling, M_alphaf, M_alphas);
                M_robinCoupling->GlobalAssemble();
            }
        }

        if(M_DDBlockPrec!=5 && M_DDBlockPrec<9)
            M_precMatrPtr.reset(new matrix_type(*M_monolithicMap));
    }
    super::evalResidual( disp, M_rhsFull, res, false );
}


void fullMonolithic::setupBlockPrec(vector_type& rhs)
{

    boost::shared_ptr<IfpackComposedPrec>  ifpackCompPrec;

    if(M_dataFluid->useShapeDerivatives())
    {
        if(M_DDBlockPrec==6)
        {
            M_precMatrPtr.reset(new matrix_type(*M_monolithicMatrix));
        }
        matrix_ptrtype monolithicMatrix(new matrix_type(*M_monolithicMap));
        *monolithicMatrix+=*M_monolithicMatrix;
        this->shapeDerivatives(monolithicMatrix,*M_uk /*subX*/, M_domainVelImplicit, M_convectiveTermDer);
        monolithicMatrix->GlobalAssemble();
        M_monolithicMatrix=monolithicMatrix;

        //meshMotion().updateDispDiff();//updates M_dispDiff to xk-xn-1
    }
    switch(M_DDBlockPrec)
    {
    case 1:
        {
            M_meshMotion->setMatrix(M_precMatrPtr);
        }
        break;
    case 2:
        {
            super::applyPreconditioner(M_robinCoupling, rhs);
        }
        break;
    case 5:
        {
            super::applyPreconditioner(M_robinCoupling, rhs);
            M_monolithicMatrix->GlobalAssemble();
        }
        break;
    case 6:
        {
            super::applyPreconditioner(M_robinCoupling, rhs);
            super::applyPreconditioner(M_robinCoupling, M_precMatrPtr);
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
            addDiagonalEntries(1., M_fluidBlock, M_mmFESpace->map(), mapWithoutMesh().getMap(Unique)->NumGlobalElements());
            couplingMatrix(M_fluidBlock, 15);

            M_BCh_Robin->setOffset(M_offset);
            M_BCh_flux->setOffset(M_offset-M_fluxes);
            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_flux, M_uFESpace->feBd(), 1., dataSolid().getTime() );
            bcManageMatrix( *M_fluidBlock, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataSolid().getTime() );
            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_u, M_uFESpace->feBd(), 1., dataSolid().getTime() );
            M_fluidBlock->GlobalAssemble();
            bcManageMatrix( *M_fluidBlock, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataSolid().getTime() );
            super::applyPreconditioner(M_robinCoupling, rhs);
            super::applyPreconditioner(M_robinCoupling, M_fluidBlock);
            M_fluidBlock->GlobalAssemble();
            M_fluidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_fluidBlock));

            ifpackCompPrec =  boost::dynamic_pointer_cast< IfpackComposedPrec, prec_raw_type > (M_precPtr);

            if(ifpackCompPrec->set())
            {
                ifpackCompPrec->replace(M_fluidOper, 0);
            }
            else
            {
                M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap));
                addDiagonalEntries(1., M_solidBlockPrec, mapWithoutMesh(), 0);
                couplingMatrix(M_solidBlockPrec, 16);
                M_meshMotion->setMatrix(M_solidBlockPrec);
                M_solidBlockPrec->GlobalAssemble();
                M_solidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_solidBlockPrec));

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
            addDiagonalEntries(1., M_fluidBlock, M_mmFESpace->map(), mapWithoutMesh().getMap(Unique)->NumGlobalElements());
            couplingMatrix(M_fluidBlock, 15);
            //*M_fluidBlock+=*M_SDMatrix;
            this->shapeDerivatives(M_fluidBlock,*M_uk /*subX*/, M_domainVelImplicit, M_convectiveTermDer);

            M_BCh_Robin->setOffset(M_offset);
            bcManageMatrix( *M_fluidBlock, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataSolid().getTime() );
            M_BCh_flux->setOffset(M_offset-M_fluxes);
            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_flux, M_uFESpace->feBd(), 1., dataSolid().getTime() );
            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_u, M_uFESpace->feBd(), 1., dataSolid().getTime() );

            bcManageMatrix( *M_fluidBlock, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataSolid().getTime() );

            super::applyPreconditioner(M_robinCoupling, rhs);
            super::applyPreconditioner(M_robinCoupling, M_fluidBlock);
            M_fluidBlock->GlobalAssemble();
            M_fluidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_fluidBlock));

            ifpackCompPrec =  boost::dynamic_pointer_cast< IfpackComposedPrec, prec_raw_type > (M_precPtr);

            if(ifpackCompPrec->set())
            {
                ifpackCompPrec->replace(M_fluidOper, 1);
            }
            else
            {
                M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap));
                addDiagonalEntries(1., M_solidBlockPrec, mapWithoutMesh(), 0);
                M_meshMotion->setMatrix(M_solidBlockPrec);
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
            addDiagonalEntries(1., M_fluidBlock, M_mmFESpace->map(), mapWithoutMesh().getMap(Unique)->NumGlobalElements());
            addDiagonalEntries(1., M_fluidBlock, M_dFESpace->map(), M_offset);
            couplingMatrix(M_fluidBlock, 7);
            //*M_fluidBlock+=*M_SDMatrix;
            this->shapeDerivatives(M_fluidBlock,*M_uk /*subX*/, M_domainVelImplicit, M_convectiveTermDer);

            M_BCh_flux->setOffset(M_offset-M_fluxes);
            if ( !M_BCh_flux->bdUpdateDone() )
                M_BCh_flux->bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );//to kill

            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_flux, M_uFESpace->feBd(), 1., dataSolid().getTime() );
            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_u, M_uFESpace->feBd(), 1., dataSolid().getTime() );

            M_fluidBlock->GlobalAssemble();

            M_fluidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_fluidBlock));

            ifpackCompPrec =  boost::dynamic_pointer_cast< IfpackComposedPrec, prec_raw_type > (M_precPtr);

            if(ifpackCompPrec->set())
            {
                ifpackCompPrec->replace(M_fluidOper, 2);
            }
            else
            {
                M_meshBlock.reset(new matrix_type(*M_monolithicMap));
                addDiagonalEntries(1., M_meshBlock, mapWithoutMesh(), 0);
                M_meshMotion->setMatrix(M_meshBlock);
                couplingMatrix(M_meshBlock, 16);
                M_meshBlock->GlobalAssemble();


                M_meshOper.reset(new IfpackComposedPrec::operator_raw_type(*M_meshBlock));


                M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap, 1));
                *M_solidBlockPrec += *M_solidBlock;
                addDiagonalEntries(1., M_solidBlockPrec, M_uFESpace->map() );
                addDiagonalEntries(1., M_solidBlockPrec, M_pFESpace->map()+M_fluxes , M_uFESpace->dof().numTotalDof()*nDimensions );
                addDiagonalEntries(1., M_solidBlockPrec, M_monolithicInterfaceMap, M_solidAndFluidDim);
                addDiagonalEntries(1., M_solidBlockPrec, M_mmFESpace->map(), M_solidAndFluidDim + nDimensions*M_interface);
                M_BCh_Robin->setOffset(M_offset);
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataSolid().getTime() );
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataSolid().getTime() );
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
            addDiagonalEntries(1., M_fluidBlock, M_mmFESpace->map(), mapWithoutMesh().getMap(Unique)->NumGlobalElements());
            addDiagonalEntries(1., M_fluidBlock, M_dFESpace->map(), M_offset);
            couplingMatrix(M_fluidBlock, 7);
            //this->shapeDerivatives(M_fluidBlock,*M_uk /*subX*/);

            M_BCh_flux->setOffset(M_offset-M_fluxes);
            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_flux, M_uFESpace->feBd(), 1., dataSolid().getTime() );
            bcManageMatrix( *M_fluidBlock, *M_uFESpace->mesh(), M_uFESpace->dof(), *this->M_BCh_u, M_uFESpace->feBd(), 1., dataSolid().getTime() );
            //	      M_fluidBlock->GlobalAssemble();
            M_fluidBlock->GlobalAssemble();
            //M_fluidBlock->spy("fb");
            M_fluidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_fluidBlock));

            ifpackCompPrec =  boost::dynamic_pointer_cast< IfpackComposedPrec, prec_raw_type > (M_precPtr);

            if(ifpackCompPrec->set())
            {
                ifpackCompPrec->replace(M_fluidOper, 2);
            }
            else
            {
                M_solidBlockPrec.reset(new matrix_type(*M_monolithicMap, 1));
                *M_solidBlockPrec += *M_solidBlock;
                addDiagonalEntries(1., M_solidBlockPrec, M_uFESpace->map() );
                addDiagonalEntries(1., M_solidBlockPrec, M_pFESpace->map(), M_uFESpace->dof().numTotalDof()*nDimensions );
                addDiagonalEntries(1., M_solidBlockPrec, M_monolithicInterfaceMap, M_solidAndFluidDim);
                addDiagonalEntries(1., M_solidBlockPrec, M_mmFESpace->map(), M_solidAndFluidDim + nDimensions*M_interface);
                M_BCh_Robin->setOffset(M_offset);
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *this->M_BCh_Robin, M_dFESpace->feBd(), 1., dataSolid().getTime() );
                bcManageMatrix( *M_solidBlockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1., dataSolid().getTime() );
                M_solidBlockPrec->GlobalAssemble();

                M_solidOper.reset(new IfpackComposedPrec::operator_raw_type(*M_solidBlockPrec));



                M_meshBlock.reset(new matrix_type(*M_monolithicMap));
                addDiagonalEntries(1., M_meshBlock, mapWithoutMesh(), 0);
                M_meshMotion->setMatrix(M_meshBlock);
                //*M_meshBlock += *M_SDMatrix;
                this->shapeDerivatives(M_meshBlock,*M_uk /*subX*/, M_domainVelImplicit, M_convectiveTermDer);
                M_meshBlock->GlobalAssemble();

                M_meshOper.reset(new IfpackComposedPrec::operator_raw_type(*M_meshBlock));

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

void fullMonolithic::solveJac(vector_type       &_muk,
                              const vector_type &_res,
                              const Real       _linearRelTol)
{
    super::setMatrix(*M_monolithicMatrix);
    super::solveJac(_muk, _res, _linearRelTol);
}

void fullMonolithic::initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                                 FSIOperator::solid_type::value_type::Function const& p0,
                                 FSIOperator::solid_type::value_type::Function const& d0,
                                 FSIOperator::solid_type::value_type::Function const& w0,
                                 FSIOperator::solid_type::value_type::Function const& df0 )
{
    super::initialize(u0, p0, d0, w0, df0);

    vector_type df(M_mmFESpace->map());
    M_mmFESpace->interpolate(df0, df, M_dataSolid->getTime());

    M_un->add(df, M_solidAndFluidDim+getDimInterface());

    M_meshMotion->setDisplacement(df);
    //  M_bdf->initialize_unk(*M_un);
}

void fullMonolithic::initializeMesh(vector_ptrtype fluid_dispOld)
{
    meshMotion().initialize(*fluid_dispOld);
}

#ifdef UNDEFINED //OBSOLETE
int Epetra_FullMonolithic::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
    //    Y=X;
    vector_type x(X, M_FMOper->getCouplingVariableMap(), Unique);
    vector_type y(*M_FMOper->getCouplingVariableMap(), Unique);
    //    x.spy("x");
    //    Epetra_FEVector  dz(Y.Map());//Ax-b, that is Ax since b=0.
    vector_ptrtype rhsShapeDer(new vector_type(*M_FMOper->getCouplingVariableMap(), Unique));
    vector_ptrtype meshDeltaDisp(new vector_type(M_FMOper->mmFESpace().map(), Unique));
    vector_type subX(M_FMOper->mapWithoutMesh());
    UInt offset(M_solidAndFluidDim + M_FMOper->getDimInterface());
    y = (*M_FMOper->getMatrixPtr())*x;
    //    std::cout<<x.getEpetraVector().NumGlobalElements()<<std::endl;
    //BCManage for the shape derivatives part of the matrix
    //    vector_ptrtype p;
    //    M_FMOper->fluid().diagonalizeVec(x , M_FMOper->BCh_solid()/*, p*/);
    //End of BCManage
    meshDeltaDisp->subset(x, offset);
    //    *meshDeltaDisp += M_FMOper->meshMotion().disp();
    //  subX.subset(*M_FMOper->uk(), 0);
    //M_FMOper->shapeDerivatives(rhsShapeDer , meshDeltaDisp,*M_FMOper->uk() /*subX*/);
    *rhsShapeDer *= -1.;//e-1;
    //rhsShapeDer->spy("solBefore");
    //BCManage for the shape derivatives part of the matrix
    //M_FMOper->fluid().bcManageVec(*rhsShapeDer, M_FMOper->BCh_sd());
    //End of BCManage
    //rhsShapeDer->spy("solAfter");
    y += *rhsShapeDer;
    Y=y.getEpetraVector();
    //    Y=X;
    //    Y.Update(1., X, 1.);
    return 0;
}
#endif
/*typedef Monolithic::vector_type vector_type;
  vector_type&
  fullMonolithic::meshVel()
  {
  vector_type meshv(M_meshMotion->dispDeltaDiff());
  meshv -= *M_meshVel;
  meshv *= M_dataFluid->getTimeStep();
  return meshv;
  }
*/

/*
  void
  fullMonolithic::solveJac(vector_type&       _muk,
  const vector_type& _res,
  const double       _linearRelTol)
  {
  super::solveJac(vector_type&       _muk,
  const vector_type& _res,
  const double       _linearRelTol);

  }*/
/*void
  fullMonolithic::updateSystem(const vector_type& solution)
  {
  UInt offset(M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof() + M_dFESpace->map().getMap(Unique)->NumGlobalElements() + nDimensions*M_interface);
  vector_type meshDisp(M_mmFESpace->map(), Repeated);
  M_meshMotion->initialize(meshDisp.subset(solution, offset));
  super::updateSystem(solution);
  }*/

//  namespace
//  {
//   FSIOperator* createFM(){ return new fullMonolithic(); }
//  }
//   static bool reg = FSIFactory::instance().registerProduct( "fullMonolithic", &createFM );

}
