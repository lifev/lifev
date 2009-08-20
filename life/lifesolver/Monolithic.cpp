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

#include <life/lifesolver/Monolithic.hpp>
//#include <life/lifesolver/reducedLinFluid.hpp>

namespace LifeV
{
// Constructors
Monolithic::Monolithic():
    super(),
    M_updateEvery(0),
    M_monolithicMap(),
    M_interfaceMap(),
    M_numerationInterface(),
    M_interface(0),
    M_beta(),
    M_reusePrec(true),
    M_resetPrec(true),
    M_maxIterSolver(-1),
    M_couplingMatrix(),
    M_linearSolver(),
    M_monolithicMatrix()
{
}
// Destructor
Monolithic::~Monolithic()
{
}


void Monolithic::setup()
{
    super::setup();


            // Added here the code to build the bigMap!!
            // Note: for now works only with matching grids (and poly order) on the interface
    assert(M_fluidInterfaceMap->getMap(Unique)->NumGlobalElements() == M_solidInterfaceMap->getMap(Unique)->NumGlobalElements());

                M_interfaceMap = *M_solidInterfaceMap;
            //            M_dofInterface = dofInterfaceSolid();

            std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
            std::map<ID, ID>::const_iterator ITrow;
            M_monolithicMap.reset(new EpetraMap( M_uFESpace->map()));
            *M_monolithicMap+= M_pFESpace->map();
            *M_monolithicMap+= M_dFESpace->map();

            int numtasks;
            MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
            int* numInterfaceDof(new int[numtasks]);
            int pid=M_epetraWorldComm->MyPID();
            int numMyElements = M_interfaceMap.getMap(Unique)->NumMyElements();
            numInterfaceDof[pid]=numMyElements;
            EpetraMap subMap(*M_interfaceMap.getMap(Unique), (UInt)0, M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions);

            M_numerationInterface.reset(new vector_type(subMap,Unique));//should be an int vector instead of double
            //                    M_numerationInterfaceInt.reset(new Epetra_IntVector(*M_interfaceMap.getMap(Unique)));

            for(int j=0; j<numtasks; ++j)
                M_epetraWorldComm->Broadcast( &numInterfaceDof[j], 1, j);

            for(int j=numtasks-1; j>0 ; --j)
                {
                    numInterfaceDof[j] = numInterfaceDof[j-1];
                }
            numInterfaceDof[0]=0;
            for(int j=1; j<numtasks ; ++j)
                numInterfaceDof[j] += numInterfaceDof[j-1];

            UInt k=1;
            UInt l=0;

            M_interface = M_interfaceMap.getMap(Unique)->NumGlobalElements()/nDimensions;
            UInt solidDim=M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions;
            //                    UInt localInterface=M_interfaceMap.getMap(Unique)->NumMyElements()/nDimensions;

            for(l=0, ITrow=locDofMap.begin(); ITrow!=locDofMap.end() ; ++ITrow)
                {
                    if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/)>=0)
                        {
                            (*M_numerationInterface)[ITrow->second /*+ dim*solidDim*/ ]=l+1+ (int)(numInterfaceDof[pid]/nDimensions)/*+ dim*localInterface*/      ;
                            //                                    (*M_numerationInterfaceInt)[ITrow->second /*+ dim*solidDim*/ ]=l+1+ (int)(M_numInterfaceDof[pid]/nDimensions)/*+ dim*localInterface*/      ;
                            if((int)(*M_numerationInterface)(ITrow->second )!=floor(l+1+ numInterfaceDof[pid]/nDimensions+0.2 /*+ dim*localInterface*/) )
                                std::cout<<"ERROR! the numeration of the coupling map is not correct"<<std::endl;

                            ++l;
                        }

                }

            //          M_numerationInterface->spy("numeration");

            std::vector<int> couplingVector;
            couplingVector.reserve((int)(M_interfaceMap.getMap(Unique)->NumMyElements()/nDimensions));

            for(int dim=0; dim<nDimensions; ++dim)
                {
                    for( ITrow=locDofMap.begin(); ITrow!=locDofMap.end() ; ++ITrow)
                        {
                            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second)>=0)
                                {
                                    couplingVector.push_back((*M_numerationInterface)(ITrow->second /*+ dim * solidDim*/)+ dim * M_interface );
                                    //couplingVector.push_back((*M_numerationInterfaceInt)[ITrow->second /*+ dim * solidDim*/]+ dim * M_interface );
                                }
                        }
                }// so the map for the coupling part of the matrix is just Unique

            EpetraMap newMap = EpetraMap(-1, couplingVector.size(), &couplingVector[0], M_monolithicMap->getMap(Repeated)->IndexBase()/*1*/, *M_epetraWorldComm);
            *M_monolithicMap  += newMap;

            //            std::cout<<"map global elements : "<<M_monolithicMap.getMap(Unique)->NumGlobalElements()<<std::endl;
            //            std::cout<<"map My elements : "<<M_monolithicMap.getMap(Unique)->NumMyElements()<<std::endl;

            //            std::cout<<"newmap global elements : "<<newMap.getMap(Unique)->NumGlobalElements()<<std::endl;
            //          std::cout<<"newmap My elements : "<<newMap.getMap(Unique)->NumMyElements()<<std::endl;

            //          std::cout<<"repeated newmap global elements : "<<newMap.getMap(Repeated)->NumGlobalElements()<<std::endl;
            //          std::cout<<"repeated newmap My elements : "<<newMap.getMap(Repeated)->NumMyElements()<<std::endl;


            //the map for the interface coupling matrices should be done with respect to the coarser mesh.

            M_meshMotion.reset(new FSIOperator::meshmotion_raw_type(*M_mmFESpace,
                                                                    *M_epetraComm));

            M_fluid.reset(new FSIOperator::fluid_raw_type(dataFluid(),
                                                         *M_uFESpace,
                                                         *M_pFESpace,
                                                         *M_epetraComm,
                                                         *M_monolithicMap));

//             if (isLinearFluid())// to be implemented
//                 M_fluidLin.reset(new FSIOperator::fluidlin_raw_type(dataFluid(),
//                                                                    *M_uFESpace,
//                                                                    *M_pFESpace,
//                                                                    *M_epetraComm));

            vector_type u0(*M_monolithicMap);
            M_bdf.reset(new BdfT<vector_type>(M_dataFluid->getBDF_order()));
            M_bdf->initialize_unk(u0);

            M_rhs.reset(new vector_type(*this->M_monolithicMap));
            M_un.reset (new vector_type(*this->M_monolithicMap));
            M_beta.reset  (new vector_type(/*M_monolithicMap*/M_uFESpace->map()));

            UInt offset = M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof();

            M_solid.reset(new FSIOperator::solid_raw_type(dataSolid(),
                                                         *M_dFESpace,
                                                         *M_epetraComm,
                                                         *M_monolithicMap,
                                                         offset
                                                         ));

//             if (isLinearSolid())// to be implemented with the offset
//                 M_solidLin.reset(new FSIOperator::solidlin_raw_type(dataSolid(),
//                                                                    *M_dFESpace,
//                                                                    *M_epetraComm));
}//end setup

void
Monolithic::setDataFromGetPot( GetPot const& data_file )
{
    super::setDataFromGetPot(data_file);

    M_updateEvery = data_file("problem/updateEvery", 0);
    this->M_dataFluid->setSemiImplicit( data_file("problem/semiImplicit", false) );
    this->M_dataFluid->setUseShapeDerivatives( data_file("fluid/useShapeDerivatives", false) );
    M_isDiagonalBlockPrec = data_file( "interface/diagonalBlockPrec",  false );
}

void
Monolithic::buildSystem()
{
    M_couplingMatrix.reset(new matrix_type(*M_monolithicMap/*, this->M_fluid->getMeanNumEntries()*/));
    couplingMatrix(M_couplingMatrix);
    M_couplingMatrix->GlobalAssemble();

    M_solid->buildSystem();
    M_solid->rescaleMatrices();
    updateCoupling(*M_couplingMatrix);
}

void
Monolithic::couplingMatrix(matrix_ptrtype & bigMatrix) // not working with non-matching grids
{
    UInt offset = M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof();
    UInt solidAndFluidDim(offset + M_dFESpace->dof().numTotalDof()*nDimensions);
    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;

    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
                    bigMatrix->set_mat_inc( offset + ITrow->second-1 + dim* M_dFESpace->dof().numTotalDof(),(int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + solidAndFluidDim, (1.) );//right low
                    bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + solidAndFluidDim, (offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), (-1.*M_solid->rescaleFactor()/*/M_dataFluid->getTimeStep()*/));//low right
                    bigMatrix->set_mat_inc( ITrow->first-1 + dim* M_uFESpace->dof().numTotalDof(), (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + solidAndFluidDim, (-1.) );//right up
                    bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + solidAndFluidDim, (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), 1.);//low left
                }
        }
    }
}


void
Monolithic::couplingRhs(vector_ptrtype rhs, vector_ptrtype un) // not working with non-matching grids
{
    UInt offset = M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof();
    UInt solidAndFluidDim(offset + M_dFESpace->dof().numTotalDof()*nDimensions);
    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    //    UInt solidDim=M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions;

    vector_type lambda(M_interfaceMap, Unique);
    this->monolithicToInterface(lambda, *un);

    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {

                    if(rhs.get() != 0)
                        (*rhs)[  (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] + dim*M_interface +solidAndFluidDim ] = -lambda( ITrow->second + dim*M_dFESpace->dof().numTotalDof()/*M_dFESpace->dof().numTotalDof()*/ )*(M_solid->rescaleFactor())/*/M_dataFluid->getTimeStep()*/;
                    //                    std::cout<<(int)(*this->M_numerationInterface)[ITrow->second/*+dim*solidDim*/] - 1 + dim*M_interface<<std::endl;

                }
        }
    }
}


void
Monolithic::addDiagonalEntries(Real& entry, matrix_ptrtype & bigMatrix)
{
    UInt solidAndFluidDim(M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof() + M_dFESpace->dof().numTotalDof()*nDimensions);
    UInt N=nDimensions*M_interface;
    for(UInt i=0; i<N ; ++i )
        bigMatrix->set_mat_inc( i + solidAndFluidDim , i + solidAndFluidDim, entry);
}

void
Monolithic::updateSystem(const vector_type& displacement)
{

    vector_type solution(*this->M_monolithicMap);
    monolithicToFluid(displacement, solution);
    //*M_rhs *= 0.;//useless? (done in evalResidual)

    this->M_bdf->shift_right(solution);

    M_meshMotion->updateSystem();

    //    transferMeshMotionOnFluid(M_meshMotion->disp(),
    //                                     *this->M_dispFluidMeshOld);
    *this->M_un                = displacement;
    this->fluid().updateUn(*this->M_un);
    *M_rhs*=0;

//    this->setDispSolid();
}

void
Monolithic::monolithicToFluid(const vector_type& disp, vector_type& dispFluid)
{
    if(disp.getMaptype()== Repeated)
        {
            vector_type dispUnique(disp, Unique);
            monolithicToFluid(dispUnique, dispFluid);
            dispFluid = dispUnique;
            return;
        }

            for(UInt k=1 ; k <= nDimensions*M_uFESpace->dof().numTotalDof() ; ++k)// since map iterators start from one, operator [] subtracts one!
                {
                    if(M_monolithicMap->getMap(Unique)->LID(k) >=0)
                        {
                            dispFluid[k] = disp(k);
                        }
                }

}


void
Monolithic::monolithicToInterface(vector_type& lambdaSolid, const vector_type& disp)
{
       if (disp.getMaptype() == Repeated)
        {
            vector_type const  dispUnique(disp, Unique);
            monolithicToInterface(lambdaSolid, dispUnique);
            return;
        }
       if (lambdaSolid.getMaptype() == Repeated)
        {
            vector_type  lambdaSolidUn(lambdaSolid.getMap(), Unique);
            monolithicToInterface( lambdaSolidUn, disp);
            lambdaSolid = lambdaSolidUn;
            return;
            }
       /* UInt MyOffset(M_uFESpace->map().getMap(Unique)->NumMyElements()+M_pFESpace->map().getMap(Unique)->NumMyElements());
       vector_type subDisp(this->M_dFESpace->map(), Unique);
       subDisp.mySubset(disp, MyOffset);
       lambdaSolid=subDisp;*/

       UInt offset(M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof());
       EpetraMap subMap(*disp.getMap().getMap(Unique), offset,disp.getMap().getMap(Unique)->NumGlobalElements() );
       vector_type subDisp(subMap, Unique);
       subDisp.subset(disp, offset);
       lambdaSolid=subDisp;

}


void
Monolithic::monolithicToSolid(const vector_type& disp, vector_type& dispSolid )// to use the subset method!
{
    if(disp.getMaptype()== Repeated)
        {
            vector_type dispUnique(disp, Unique);
            monolithicToSolid(dispUnique, dispSolid);
            return;
        }
    if (dispSolid.getMaptype() == Repeated)
        {
            vector_type  dispSolidUn(dispSolid.getMap(), Unique);
            monolithicToSolid( disp,  dispSolidUn);
            dispSolid=dispSolidUn;
            return;
        }

    dispSolid *= 0.;
    UInt offset(M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof());
    UInt k = 0;
    UInt N = M_dFESpace->dof().numTotalDof();

    for( UInt dim = 0; dim<nDimensions; ++dim)
        {
            for( k=1 ; k <= N ; ++k)// since map iterators start from one, operator [] subtracts one!
                {
                    if(M_monolithicMap->getMap(Unique)->LID(k+N*dim+offset) >= 0 )
                        {
                            dispSolid[k + N*dim + offset] = disp(k + N*dim + offset);
                        }
                }
        }
}

void Monolithic::setDispSolid(const vector_type &sol)
{
    vector_type disp(*M_monolithicMap);
    monolithicToSolid(sol, disp);
    this->M_solid->setDisp(disp);
}

void
Monolithic::evalResidual( vector_type&       res,
                          const vector_type& disp,
                          const UInt          iter )
{

            //            eval(disp, iter);

    Chrono chronoFluid, chronoSolid, chronoInterface;

    setDispSolid(disp);
    M_nbEval++ ;

    if(iter == 0)
        {
            M_nbEval = 0; // new time step
            this->M_solid->resetPrec();
        }

    if(iter==0 || !this->M_dataFluid->isSemiImplicit())
        {

                vector_type lambdaFluid(this->M_interfaceMap, Unique);

                monolithicToInterface(lambdaFluid, disp);

                lambdaFluid *= M_dataFluid->getTimeStep();

                this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface

                //                M_epetraWorldComm->Barrier();
                chronoFluid.start();

                M_meshMotion->iterate();
                M_meshMotion->updateDispDiff();
                //                this->transferMeshMotionOnFluid(M_meshMotion->disp(),
                //                                        this->veloFluidMesh());

                //                this->veloFluidMesh()    -= dispFluidMeshOld();
                //                this->veloFluidMesh()    *= 1./(M_dataFluid->getTimeStep());

                // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
                // the value of the displacement for some points
                                *this->M_beta *= 0.;
                                //*this->M_rhs *= 0.;

                vector_type const meshDispDiff( M_meshMotion->dispDiff(), Repeated );

                this->moveMesh(meshDispDiff);

                this->interpolateVelocity(meshDispDiff, *this->M_beta);

                *this->M_beta *= -1./M_dataFluid->getTimeStep();
                vector_ptrtype fluid(new vector_type(this->M_uFESpace->map()));
                fluid->subset(*M_un, 0);
                *this->M_beta += *fluid/*M_un*/;

                double alpha = 1./M_dataFluid->getTimeStep();

                M_monolithicMatrix.reset(new matrix_type(*M_monolithicMap/*, this->M_fluid->getMeanNumEntries()*/));

                M_fluid->updateSystem(alpha,* this->M_beta, *this->M_rhs );//here it assembles the fluid matrices
                // coupling matrices assebling
                M_fluid->updateStab( *M_monolithicMatrix);//applies the stabilization terms
                M_fluid->getFluidMatrix( *M_monolithicMatrix);
                updateMatrix( *M_monolithicMatrix);

                if(iter==0)
                    {
                        M_nbEval = 0; // new time step
                        *this->M_rhs   += M_fluid->matrMass()*M_bdf->time_der( M_dataFluid->getTimeStep() );
                        couplingRhs( this->M_rhs, M_un);
                        this->M_solid->updateStuff();
                        updateSolidSystem(this->M_rhs);
                    }


                //*monolithicMatrix += *this->M_fluid->getBigMatrixPtr();


                //M_solid->updateMatrix( *monolithicMatrix);

                //if(iter==0)


                    //this->setDispSolid();

                //this->updateSolidSystem(this->M_rhs);

                this->evalResidual( *M_BCh_u, *M_BCh_d, disp, M_rhs, res);
                //this->evalResidual( *M_BCh_u, *M_BCh_d, disp, res/*,  *this->M_fluid->getBigMatrixPtr() */);
            }

        vector_ptrtype residual(new vector_type(*M_monolithicMap, Unique));
        this->evalResidual( disp, M_rhs, *residual );
        res = *residual;

}

void Monolithic::
evalResidual( const vector_type& sol,const vector_ptrtype& rhs, vector_type& res, bool diagonalScaling)
{
    if(diagonalScaling)
        diagonalScale(*rhs, M_monolithicMatrix);
    res = *M_monolithicMatrix*sol;//solid->getMatrixPtr()*sol;
    res -= *rhs;
    //    rhs->spy("rhsAfter");
    // Ax-b
}

void Monolithic::
evalResidual( fluid_bchandler_raw_type& bchFluid, solid_bchandler_raw_type& bchSolid, const vector_type& sol, const vector_ptrtype& rhs, vector_type& res, bool diagonalScaling, matrix_ptrtype preconditioner)
{
    //M_solid->getMatrixPtr()->GlobalAssemble();
    M_monolithicMatrix->GlobalAssemble();
    matrix_ptrtype tmpMatPtr(new matrix_type(*M_monolithicMatrix));
    tmpMatPtr->GlobalAssemble();
    //preconditioner->GlobalAssemble();
    //M_solid->getMatrixPtr().reset(new matrix_type(M_solid->getMap()));
    M_monolithicMatrix.reset(new matrix_type(M_solid->getMap()));
    int err = err = EpetraExt::MatrixMatrix::Multiply( preconditioner->getEpetraMatrix(), false, tmpMatPtr->getEpetraMatrix(), false, M_monolithicMatrix->getEpetraMatrix());
    *rhs = (*preconditioner)*(*rhs);
    evalResidual( bchFluid, bchSolid, sol, rhs, res, diagonalScaling);
}

void Monolithic::
evalResidual( fluid_bchandler_raw_type& bchFluid, solid_bchandler_raw_type& bchSolid, const vector_type& sol, const vector_ptrtype& rhs, vector_type& res, bool diagonalScaling)
{
    vector_type rhsFullSolid(*rhs, Unique); // ignoring non-local entries, Otherwise they are summed up lately
    //    applyBoundaryConditions( *M_solid->getMatrixPtr(), rhsFullSolid, bchSolid, M_solid->offset() );
        if(M_solid->offset())
            bchSolid.setOffset(M_solid->offset());
        if ( !bchSolid.bdUpdateDone() )
            bchSolid.bdUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );

        bcManage( *M_monolithicMatrix, rhsFullSolid, *M_dFESpace->mesh(), M_dFESpace->dof(), bchSolid, M_dFESpace->feBd(), 1.,
                  dataSolid().getTime() );

    // matrix should be GlobalAssembled by  bcManage

    if ( !bchFluid.bdUpdateDone() )
        bchFluid.bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );

    vector_type rhsFull(rhsFullSolid, Unique);
    bcManage( *M_monolithicMatrix, rhsFull, *M_uFESpace->mesh(), M_uFESpace->dof(), bchFluid, M_uFESpace->feBd(), 1.,
              dataSolid().getTime() );

    M_solid->getDisplayer().leaderPrint("rhs norm = ", rhs->NormInf() );

    *rhs = rhsFull;

    evalResidual(sol,rhs, res, diagonalScaling);
}




void  Monolithic::solveJac(vector_type         &_step,
                           const vector_type   &_res,
                           const double         /*_linearRelTol*/)
{

    std::cout << "solveJac: NormInf res " << _res.NormInf() << std::endl;
    _step *= 0.;

    std::cout << "Solving Jacobian system... " << std::endl;
    //    M_epetraWorldComm->Barrier();

    matrix_ptrtype bigPrecPtr;

    if(this->M_fluid->getIsDiagonalBlockPrec()==true /*M_diagonalBlockPreconditioner*/)
        {
            bigPrecPtr.reset(new matrix_type(*M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
            this->M_fluid->setBlockPreconditioner(bigPrecPtr);
                //            *bigPrecPtr += *this->M_fluid->getBlockPrecPtr();
            this->M_fluid->updateStab(*bigPrecPtr);
            Real entry(1.0);
            this->setBlockPreconditioner(bigPrecPtr);
            addDiagonalEntries(entry,bigPrecPtr);
            if(!M_isDiagonalBlockPrec)
                couplingMatrix(bigPrecPtr);
            bigPrecPtr->GlobalAssemble();
        }
    else
        {
            if(M_isDiagonalBlockPrec)
                {
                    bigPrecPtr.reset(new matrix_type(*M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
                    this->M_fluid->updateStab( *bigPrecPtr);//applies the stabilization terms
                    this->M_fluid->getFluidMatrix( *bigPrecPtr);
                    this->setBlockPreconditioner(bigPrecPtr);
                    Real entry(1.0);
                    addDiagonalEntries(entry,bigPrecPtr);
                    bigPrecPtr->GlobalAssemble();
                }
            else
                {
                    bigPrecPtr.reset(new matrix_type(*M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
                    this->setFullPreconditioner(bigPrecPtr);
                    Real entry(1.0);
                    addDiagonalEntries(entry,bigPrecPtr);
                    bigPrecPtr->GlobalAssemble();
                }
        }
    this->iterateMonolithic(*(const_cast<vector_type*>(&_res))/*just to avoid the const type*/, _step, bigPrecPtr, M_linearSolver);

    //setDispSolid();// to be done just if not semiImplicit (done once yet in updateSystem)
    //this->M_solid->updateStuff();
    std::cout << "done." << std::endl;

}

void
Monolithic::iterateMesh(const vector_type& disp)
{
                vector_type lambdaFluid(this->M_interfaceMap, Unique);

                monolithicToInterface(lambdaFluid, disp);

                lambdaFluid *= M_dataFluid->getTimeStep();

                this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface

                //                M_epetraWorldComm->Barrier();

                M_meshMotion->iterate();

                //                this->transferMeshMotionOnFluid(M_meshMotion->disp(),
                //                                        this->veloFluidMesh());

                //                this->veloFluidMesh()    -= dispFluidMeshOld();
                //                this->veloFluidMesh()    *= 1./(M_dataFluid->getTimeStep());

                // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
                // the value of the displacement for some points
                //                vector_type const meshDispDiff( M_meshMotion->dispDiff(), Repeated );

                //                this->moveMesh(meshDispDiff);
}

void Monolithic::
updateCoupling(matrix_type couplingMatrix)
{
    M_couplingMatrix.reset(new matrix_type(M_solid->getMap()));
    *M_couplingMatrix += *M_solid->getMassStiff();
    *M_couplingMatrix += couplingMatrix;
    M_couplingMatrix->GlobalAssemble();
    //M_massStiff.reset(new matrix_type(*tmp));
}

void Monolithic::
applyPreconditioner( matrix_ptrtype robinCoupling)
{
    matrix_type tmpMatrix(*M_monolithicMatrix);
    tmpMatrix.GlobalAssemble();
    M_monolithicMatrix.reset(new matrix_type(M_solid->getMap()));
    EpetraExt::MatrixMatrix::Multiply( robinCoupling->getEpetraMatrix(), false, tmpMatrix.getEpetraMatrix(), false, M_monolithicMatrix->getEpetraMatrix());
    *M_rhs=*robinCoupling*(*M_rhs);
}


void Monolithic::
updateMatrix(matrix_type & bigMatrixStokes)
{
    *M_monolithicMatrix += *M_couplingMatrix;
    M_monolithicMatrix->GlobalAssemble();
}

void Monolithic::
updateSolidSystem( vector_ptrtype & rhsFluidCoupling )
{
    M_solid->updateSystem();
    *rhsFluidCoupling += *M_solid->rhsWithoutBC();
}
void
Monolithic::setUpSystem( GetPot const& data_file )
{
        M_fluid->setUp(data_file);
        M_meshMotion->setUp(data_file);
        setUp(data_file);
}

void Monolithic::
setUp( const GetPot& dataFile )
{
    M_solid->getDisplayer().leaderPrint("\n S-  Displacement unknowns: ",  M_dFESpace->dof().numTotalDof() );
    M_solid->getDisplayer().leaderPrint(" S-  Computing mass and linear strain matrices ... \n");
    M_linearSolver.reset(new SolverTrilinos(*M_epetraComm));
    M_reusePrec     = dataFile( "solid/prec/reuse", true);
    //M_maxIterSolver = dataFile( "solid/solver/max_iter", -1);

    M_linearSolver->setDataFromGetPot( dataFile, "solid/solver" );
    M_linearSolver->setUpPrec(dataFile, "solid/prec");
    M_reusePrec     = dataFile( "solid/prec/reuse", true);
    //    M_maxIterForReuse = data_file( "solid/solver/max_iter_reuse", M_maxIterSolver*8/10);
    M_maxIterSolver = dataFile( "solid/solver/max_iter", -1);
}

void Monolithic::
diagonalScale(vector_type& rhs, matrix_ptrtype matrFull)
{
    Epetra_Vector diagonal(*rhs.getMap().getMap(Unique));
    //M_matrFull->getEpetraMatrix().InvRowSums(diagonal);
    //M_matrFull->getEpetraMatrix().InvRowMaxs(diagonal);
    //M_matrFull->getEpetraMatrix().InvColSums(diagonal);
    matrFull->getEpetraMatrix().InvColMaxs(diagonal);
    matrFull->getEpetraMatrix().LeftScale(diagonal);
    rhs.getEpetraVector().Multiply(1, rhs.getEpetraVector(), diagonal,0);
}


//void Monolithic::solidInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct)
void Monolithic::solidInit(const std::string& dOrder)
{   // Monolitic: In the beginning I need a non-partitioned mesh. later we will do the partitioning
    M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(M_dataSolid->mesh(),
														dOrder,
                                                       //*refFE_struct,
                                                       //*qR_struct,
                                                       //*bdQr_struct,
                                                       3,
                                                       *M_epetraComm));
}

//void Monolithic::variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct)
void Monolithic::variablesInit(const std::string& dOrder)
{
    //EpetraMap interfaceMap(*M_solidInterfaceMap);
    //   M_solidMeshPart.reset( new  partitionMesh< FSIOperator::mesh_type > (*M_dataSolid->mesh(), *M_epetraComm, M_solidInterfaceMap->getMap(Unique).get(), M_solidInterfaceMap->getMap(Repeated).get()));

    //             M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(*M_solidMeshPart,
    //                                                                *refFE_struct,
    //                                                                *qR_struct,
    //                                                                *bdQr_struct,
    //                                                                3,
    //                                                                *M_epetraComm));

    M_solidMeshPart.reset( new  partitionMesh< mesh_type > (*M_dataSolid->mesh(), *M_epetraComm));

    M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(*M_solidMeshPart,
														dOrder,
                                                       //*refFE_struct,
                                                       //*qR_struct,
                                                       //*bdQr_struct,
                                                       3,
                                                       *M_epetraComm));
    // INITIALIZATION OF THE VARIABLES
    M_lambdaFluid.reset(new vector_type(*M_fluidInterfaceMap, Unique) );
    M_lambdaFluidRepeated.reset(new vector_type(*M_fluidInterfaceMap, Repeated) );

}


namespace
{
FSIOperator* createM(){ return new Monolithic(); }
}
static bool reg = FSIFactory::instance().registerProduct( "monolithic", &createM );

}

#endif
