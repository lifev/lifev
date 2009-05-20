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
    M_linearSolver()
    //    M_meshOld(new vector_type(M_mmFESpace->map()))
    //    M_epetraOper()
{}

void
fullMonolithic::setup()
{
    super::setup();
    UInt offset = M_monolithicMap->getMap(Unique)->NumGlobalElements();
    M_mapWithoutMesh.reset(new EpetraMap(*M_monolithicMap));
    M_BCh_mesh->setOffset(offset);
    *this->M_monolithicMap += this->M_mmFESpace->map();
    /*
    if(M_dataFluid->useShapeDerivatives())
        {
            M_epetraOper.reset( new Epetra_FullMonolithic(this));
            M_solid->setOperator(*M_epetraOper);
            }*/
    vector_type u0(*M_monolithicMap);
    M_bdf.reset(new BdfT<vector_type>(M_dataFluid->order_bdf()));
    M_bdf->initialize_unk(u0);
    this->M_rhs.reset(new vector_type(*this->M_monolithicMap));
    this->M_rhsFull.reset(new vector_type(*this->M_monolithicMap));
    M_un.reset (new vector_type(*this->M_monolithicMap));
    M_unOld.reset(new vector_type(*this->M_monolithicMap));//*********************
    //    UInt offset = M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof();
    M_meshMotion.reset(new FSIOperator::meshmotion_raw_type(*M_mmFESpace,
                                                            *M_epetraComm,
                                                            *M_monolithicMap,
                                                            offset));
    M_fluid.reset(new FSIOperator::fluid_raw_type(dataFluid(),
                                                  *M_uFESpace,
                                                  *M_pFESpace,
                                                  *M_epetraComm,
                                                  *M_monolithicMap));
    offset = M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof();
    //M_solid.reset(dynamic_cast<solidInterface*>(&(*M_solid)));
//     *(dynamic_cast<solidInterface*>(&(*M_solid)))= solid_raw_type(dataSolid(),
//                                                   *M_dFESpace,
//                                                   *M_epetraComm,
//                                                   *M_monolithicMap,
//                                                   offset,
//                                                   M_uFESpace
//                                                   );
     M_solid.reset(new solid_raw_type(dataSolid(),
                                                   *M_dFESpace,
                                                   *M_epetraComm,
                                                   *M_monolithicMap,
                                                   offset
                                                   ));

}

//fullMonolithic::fullMonolithic():
//    super()
//{
//}


void
fullMonolithic::updateSystem(const vector_type& solution)
{
    //    M_meshVel.reset(new vector_type(M_meshMotion->dispDeltaDiff()));
    *M_unOld=*M_un;
    super::updateSystem(solution);
    M_uk.reset(new vector_type(solution));

//     UInt offset(M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof() + M_dFESpace->dof().numTotalDof()*nDimensions + nDimensions*M_interface);

//     vector_ptrtype meshDispDiff(new vector_type(M_mmFESpace->map()));
//     meshDispDiff->subset(solution, offset);
//     M_meshMotion->initialize(*meshDispDiff);//M_disp is set to the total mesh disp.
}

void
fullMonolithic::couplingMatrix(matrix_ptrtype         &fullMMatrix, bool solidCoupling)
{
    super::couplingMatrix(fullMMatrix, solidCoupling);
    UInt offset(M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof());
    UInt solidFluidInterface(offset + M_dFESpace->dof().numTotalDof()*nDimensions + nDimensions*M_interface);
    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
                    fullMMatrix->set_mat_inc(solidFluidInterface + ITrow->first + dim*M_mmFESpace->dof().numTotalDof() - 1, offset + ITrow->second-1 + dim* M_dFESpace->dof().numTotalDof(), (-1.0)*M_dataFluid->timestep()*M_solid->rescaleFactor()/**1.e-2*//*scaling of the solid matrix*/ );
                }
        }
    }
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
                                const int          iter )
{
    /**
       Here I calculate the solution of the monolithic given the harmonic extension sol, then I solve the harmonic extension problem. The residual is the difference between the new found domain displacement and the old one.
    */
    setDispSolid(disp);// to be done just if not semiImplicit (done once yet in updateSystem)

    if(iter == 0)
        {
            M_nbEval = 0; // new time step
            //            this->M_solid->resetPrec();
        }
    else
        {
            this->M_solid->updateStuff();
        }
    if((iter==0)|| !this->M_dataFluid->isSemiImplicit())
        {
            this->M_beta.reset(new vector_type(M_uFESpace->map()));
//             if(iter == 0)
//                 {
//                 this->M_rhs.reset(new vector_type(*M_monolithicMap)) ;
//                 this->M_rhsFull.reset(new vector_type(*this->M_monolithicMap));
//                }

            //            M_dispOld.reset(new vector_type(M_meshMotion->dispOld()));
            UInt offset(M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof() + M_dFESpace->dof().numTotalDof()*nDimensions + nDimensions*M_interface);

            vector_ptrtype meshDispDiff(new vector_type(M_mmFESpace->map()));

            meshDispDiff->subset(disp, offset); //if the conv. term is to be condidered implicitly
            //meshDispDiff->subset(*M_uk, offset); //if the mesh motion is at the previous time step in the convective term
            //meshDispDiff->subset(*M_un, offset); //if we linearize in a semi-implicit way

            M_meshMotion->initialize(*meshDispDiff);//M_disp is set to the total mesh disp.
            double alpha = 1/M_dataFluid->timestep();


            //meshDispDiff.reset(new vector_type(M_meshMotion->disp(), Repeated));//just to repeat it

            *meshDispDiff -= M_meshMotion->dispOld();

            vector_type mmRep(*meshDispDiff, Repeated);// just to repeat dispDiff. No way witout copying?
            this->moveMesh(mmRep);
            //mmRep *= alpha;

            //*************ADDED****************
            mmRep.subset(*M_un, offset); //if the mesh motion is at the previous time step in the convective term
            mmRep -= M_meshMotion->dispOld();
            mmRep *= alpha;
            //****************END OF ADDED*********************

            this->interpolateVelocity(mmRep, *this->M_beta);
            //            *this->M_beta *= -alpha; //HE solution scaled!
            vector_ptrtype fluid(new vector_type(this->M_uFESpace->map()));
            fluid->subset(*M_un/**M_unOld*/, 0);
            *this->M_beta += *fluid/*M_un*/;
            //          if(firstIter)
            M_fluid->recomputeMatrix(true);
            M_fluid->updateSystem(alpha, *this->M_beta, *this->M_rhs );//here it assembles the fluid matrices
            if(iter==0)
                *this->M_rhs               += this->M_fluid->matrMass()*this->M_bdf->time_der( M_dataFluid->timestep() );
            M_monolithicMatrix.reset(new matrix_type(*M_monolithicMap/*, this->M_fluid->getMeanNumEntries()*/));
            M_fluid->getFluidMatrix( *M_monolithicMatrix);
            M_fluid->updateStab( *M_monolithicMatrix);//applies the stabilization terms

            /*            if(firstIter)
                          {*/
            if(iter == 0)
                super::couplingRhs( this->M_rhs, this->M_un);
                    /*
                    M_solid->buildSystem();
                    M_solid->rescaleMatrices();
                    M_solid->updateCoupling(*M_couplingMatrix);
                    firstIter=false;

                    M_meshMotion->computeMatrix();
                    }*/
            //          M_meshMotion->rescaleMatrix(dt);
            this->updateMatrix(*M_monolithicMatrix);
            //            this->setDispSolid(disp);
            if(iter == 0)
                {
                    super::updateSolidSystem(this->M_rhs);
                }
            M_uk.reset(new vector_type(disp));
            this->M_rhsFull.reset(new vector_type(*this->M_rhs));
            M_meshMotion->applyBoundaryConditions(*M_rhsFull);
            M_meshMotion->setMatrix(M_monolithicMatrix);

            super::evalResidual( *M_BCh_u, *M_BCh_d, disp, this->M_rhsFull, res, false);
            //////////////////Part under development////////////////////
//            if(M_DDBlockPrec>=2 && M_DDBlockPrec<5)
//             {
//                 if(!M_robinCoupling.get())
//                     {
//                         M_robinCoupling.reset( new matrix_type(*M_monolithicMap));
//                         super::robinCoupling(M_robinCoupling, M_alphaf, M_alphas);
//                         M_robinCoupling->GlobalAssemble();
//                     }
//                 super::applyPreconditioner(M_robinCoupling);
//                 //this->M_solid->evalResidual( disp, res, false);
//            }
             ///////////////// end///////////////////

            if(M_DDBlockPrec<5)
                M_bigPrecPtr.reset(new matrix_type(*M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
            //M_monolithicMatrix->spy("monolithicMatrix");
        }
    //else
    //{
    super::evalResidual( disp, M_rhsFull, res, false );
            //}

}

void fullMonolithic::solveJac(vector_type       &_muk,
              const vector_type &_res,
              const double       _linearRelTol)
{
    //    M_bigPrecPtr.reset(new matrix_type(M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
    /*    this->M_solid->setFullPreconditioner(M_bigPrecPtr);
    Real entry(1.e-3);
    addDiagonalEntries(entry,M_bigPrecPtr);
    M_bigPrecPtr->GlobalAssemble();
    M_bigPrecPtr->spy("bigPrecPtr");*/
    vector_ptrtype rhs(new vector_type(_res));

    if(M_dataFluid->useShapeDerivatives())
        {
            //std::cout<<"this->solidInterfaceMap()->getMap(Repeated)"<<this->solidInterfaceMap()->getMap(Repeated)<<std::endl;
            M_epetraOper.reset( new Epetra_FullMonolithic(*this/*, *this->monolithicMap()->getMap(Unique)*/));
            super::setOperator(*M_epetraOper);
            //            vector_type meshDeltaDisp(M_mmFESpace->map());
            //            UInt offset(uFESpace().dof().numTotalDof()*nDimensions +pFESpace().dof().numTotalDof() +dFESpace().map().getMap(Unique)->NumGlobalElements() + dimInterface());
            //            meshDeltaDisp.subset(*M_uk, offset);
            meshMotion().updateDispDiff();//updates M_dispDiff to xk-xn
            //BC shape Derivatives part
            /*            vector_type bcVector(*monolithicMap(), Unique);//zero vector
            vector_ptrtype rhsShapeDer(new vector_type(*monolithicMap(), Unique));
            vector_ptrtype meshDeltaDisp(new vector_type(mmFESpace().map(), Unique));
            vector_type subX(mapWithoutMesh());
            UInt offset(mapWithoutMesh().getMap(Unique)->NumGlobalElements());
            fluid().bcManageVec(bcVector, BCh_harmonicExtension());//bcVector nonzero just on the boundary
            std::cout<<"bcVector==>"<<bcVector.NormInf()<<std::endl;
            meshDeltaDisp->subset(bcVector, offset);
            subX.subset(*uk(), 0);
            shapeDerivatives(rhsShapeDer , meshDeltaDisp, subX);

            *rhs += *rhsShapeDer;*/
            //            vector_ptrtype rhs(&res);
            //            fluid().diagonalizeVec(*rhsShapeDer, M_FMOper->BCh_fluid(), rhs);
            //end of BC part
        }
    switch(M_DDBlockPrec)
        {
        case 1:
            {
                M_meshMotion->setMatrix(M_bigPrecPtr);
            }
        case 2:
            {
            }
            break;
        default:
            {}
            break;
        }

    super::solveJac(_muk, *rhs, _linearRelTol);
    //    *M_uk += _muk;
    //    UInt offset(uFESpace().dof().numTotalDof()*nDimensions +pFESpace().dof().numTotalDof() +dFESpace().map().getMap(Unique)->NumGlobalElements() + dimInterface());
    //    vector_type meshIncrement(M_mmFESpace->map());
    //    meshIncrement.subset(_muk, offset);
    //    M_meshMotion->updateDeltaDisp(meshIncrement);
}

int Epetra_FullMonolithic::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
    //    Y=X;
    vector_type x(X, M_FMOper->couplingVariableMap(), Unique);
    vector_type y(*M_FMOper->couplingVariableMap(), Unique);
    //    x.spy("x");
    //    Epetra_FEVector  dz(Y.Map());//Ax-b, that is Ax since b=0.
    vector_ptrtype rhsShapeDer(new vector_type(*M_FMOper->couplingVariableMap(), Unique));
    vector_ptrtype meshDeltaDisp(new vector_type(M_FMOper->mmFESpace().map(), Unique));
    vector_type subX(M_FMOper->mapWithoutMesh());
    UInt offset(M_FMOper->uFESpace().dof().numTotalDof()*nDimensions +  M_FMOper->pFESpace().dof().numTotalDof() + M_FMOper->dFESpace().map().getMap(Unique)->NumGlobalElements() + M_FMOper->dimInterface());
    y = (*M_FMOper->getMatrixPtr())*x;
    //    std::cout<<x.getEpetraVector().NumGlobalElements()<<std::endl;
    //BCManage for the shape derivatives part of the matrix
    //    vector_ptrtype p;
    //    M_FMOper->fluid().diagonalizeVec(x , M_FMOper->BCh_solid()/*, p*/);
    //End of BCManage
    meshDeltaDisp->subset(x, offset);
    //    *meshDeltaDisp += M_FMOper->meshMotion().disp();
    //  subX.subset(*M_FMOper->uk(), 0);
    M_FMOper->shapeDerivatives(rhsShapeDer , meshDeltaDisp,*M_FMOper->uk() /*subX*/);
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
/*typedef Monolithic::vector_type vector_type;
vector_type&
fullMonolithic::meshVel()
{
    vector_type meshv(M_meshMotion->dispDeltaDiff());
    meshv -= *M_meshVel;
    meshv *= M_dataFluid->timestep();
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
