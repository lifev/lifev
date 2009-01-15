/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

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

#include <life/lifesolver/Monolithic.hpp>
//#include <life/lifesolver/reducedLinFluid.hpp>

namespace LifeV
{
// Constructors
Monolithic::Monolithic():
    super(),
    M_updateEvery(0),
    M_monolithicMap(),
    firstIter(true),
    M_interfaceMap(),
    M_numerationInterface(),
    M_interface(0),
    M_beta()

    //    M_semiImplicit(1)
    //    M_numerationInterfaceInt()
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
            M_monolithicMap = M_uFESpace->map();
            M_monolithicMap+= M_pFESpace->map();
            M_monolithicMap+= M_dFESpace->map();

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

            EpetraMap newMap = EpetraMap(-1, couplingVector.size(), &couplingVector[0], M_monolithicMap.getMap(Repeated)->IndexBase()/*1*/, *M_epetraWorldComm);
            M_monolithicMap  += newMap;

            //            std::cout<<"map global elements : "<<M_monolithicMap.getMap(Unique)->NumGlobalElements()<<std::endl;
            //            std::cout<<"map My elements : "<<M_monolithicMap.getMap(Unique)->NumMyElements()<<std::endl;

            //            std::cout<<"newmap global elements : "<<newMap.getMap(Unique)->NumGlobalElements()<<std::endl;
            //          std::cout<<"newmap My elements : "<<newMap.getMap(Unique)->NumMyElements()<<std::endl;

            //          std::cout<<"repeated newmap global elements : "<<newMap.getMap(Repeated)->NumGlobalElements()<<std::endl;
            //          std::cout<<"repeated newmap My elements : "<<newMap.getMap(Repeated)->NumMyElements()<<std::endl;


            //the map for the interface coupling matrices should be done with respect to the coarser mesh.

            M_fluid.reset(new FSIOperator::fluid_raw_type(dataFluid(),
                                                         *M_uFESpace,
                                                         *M_pFESpace,
                                                         *M_epetraComm,
                                                         M_monolithicMap));

//             if (isLinearFluid())// to be implemented
//                 M_fluidLin.reset(new FSIOperator::fluidlin_raw_type(dataFluid(),
//                                                                    *M_uFESpace,
//                                                                    *M_pFESpace,
//                                                                    *M_epetraComm));

            vector_type u0(M_monolithicMap);
            M_bdf.reset(new BdfT<vector_type>(M_dataFluid->order_bdf()));
            M_bdf->initialize_unk(u0);

            M_rhs.reset(new vector_type(this->M_monolithicMap));
            M_un.reset (new vector_type(this->M_monolithicMap));
            M_beta.reset  (new vector_type(/*M_monolithicMap*/M_uFESpace->map()));

            UInt offset = M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof();

            M_solid.reset(new FSIOperator::solid_raw_type(dataSolid(),
                                                         *M_dFESpace,
                                                         *M_epetraComm,
                                                         M_monolithicMap,
                                                         offset,
                                                         M_uFESpace
                                                         ));

//             if (isLinearSolid())// to be implemented with the offset
//                 M_solidLin.reset(new FSIOperator::solidlin_raw_type(dataSolid(),
//                                                                    *M_dFESpace,
//                                                                    *M_epetraComm));
}//end setup

void
Monolithic::setDataFromGetPot( GetPot const& data_file )
{
//     M_solverTrilinos.setOptionsFromGetPot(data_file,"jacobian/aztec");

  super::setDataFromGetPot(data_file);

    //  this->M_dataFluid.reset(new data_fluid(data_file));
    //    this->M_dataSolid.reset(new data_solid(data_file));
    M_updateEvery = data_file("problem/updateEvery", 0);
    this->M_dataFluid->setSemiImplicit( data_file("problem/semiImplicit", false) );
    this->M_dataFluid->setUseShapeDerivatives( data_file("fluid/useShapeDerivatives", false) );
    M_isDiagonalBlockPrec = data_file( "interface/diagonalBlockPrec",  false );
}

void
Monolithic::buildSystem(){}
/*{
    if (this->isFluid())
    {
        this->M_fluid->buildSystem();
    }
    if (this->isSolid())
    {
        coupling(this->M_fluid->getBigMatrixPtr(), this->M_rhs);
        this->M_solid->buildSystem( *this->M_fluid->getBigMatrixPtr());
    }
    }*/

void
Monolithic::coupling(matrix_ptrtype & bigMatrix, vector_ptrtype rhs) // not working with non-matching grids
{

    UInt offset(M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof());
    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    UInt solidAndFluidDim(offset + M_dFESpace->dof().numTotalDof()*nDimensions);
    //    UInt solidDim=M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions;

    vector_type lambda(M_interfaceMap, Unique);
    this->monolithicToInterface(lambda, *M_un);

    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
                    bigMatrix->set_mat_inc( offset + ITrow->second-1 + dim* M_dFESpace->dof().numTotalDof(),(int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + solidAndFluidDim, (1.0) );
                    bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + solidAndFluidDim, (offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), (-1/*/M_dataFluid->timestep()*/));
                    bigMatrix->set_mat_inc( ITrow->first-1 + dim* M_uFESpace->dof().numTotalDof(), (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + solidAndFluidDim, (-1.0) );
                    bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + solidAndFluidDim, (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), 1.0);
                    bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + solidAndFluidDim , (int)(*this->M_numerationInterface)[ITrow->second /*+ dim*solidDim*/ ] - 1 + dim*M_interface + solidAndFluidDim, 0.0);
                    if(rhs.get() != 0)
                       (*rhs)[  (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] + dim*M_interface +solidAndFluidDim ] = -lambda( ITrow->second + dim*M_dFESpace->dof().numTotalDof()/*M_dFESpace->dof().numTotalDof()*/ )/*/M_dataFluid->timestep()*/;
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

    vector_type solution(this->M_monolithicMap);
    monolithicToFluid(displacement, solution);
    *M_rhs *= 0.;//useless? (done in evalResidual)

    this->M_bdf->shift_right(solution);

    M_meshMotion->updateSystem();

    //    transferMeshMotionOnFluid(M_meshMotion->disp(),
    //                                     *this->M_dispFluidMeshOld);
    *this->M_un                = displacement;

//    this->setDispSolid();
}

void
Monolithic::monolithicToFluid(const vector_type& disp, vector_type& dispFluid)
{
    if(disp.getMaptype()== Repeated)
        {
            vector_type dispUnique(disp, Unique);
            monolithicToFluid(dispUnique, dispFluid);
            return;
        }
       if (dispFluid.getMaptype() == Unique)
           /*        {
            vector_type  dispFluidRep(dispFluid.getMap(), Repeated);
            monolithicToFluid(disp, dispFluidRep);
            dispFluid = dispFluidRep;
            return;
            }*/

            for(UInt k=1 ; k <= nDimensions*M_uFESpace->dof().numTotalDof() ; ++k)// since map iterators start from one, operator [] subtracts one!
                {
                    if(M_monolithicMap.getMap(Unique)->LID(k) >=0)
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
                    if(M_monolithicMap.getMap(Unique)->LID(k+N*dim+offset) >= 0 )
                        {
                            dispSolid[k + N*dim + offset] = disp(k + N*dim + offset);
                        }
                }
        }
}

void Monolithic::setDispSolid()
{
    vector_type disp(this->M_monolithicMap);
    monolithicToSolid(this->M_solid->disp(), disp);
    this->M_solid->setDispSolid(disp);
}

void
Monolithic::evalResidual( vector_type&       res,
                          const vector_type& disp,
                          const int          iter )
{

            //            eval(disp, iter);

    Chrono chronoFluid, chronoSolid, chronoInterface;

    if(iter == 0)
        {
            M_nbEval = 0; // new time step
            this->M_solid->resetPrec();
        }

    M_nbEval++ ;

        if((iter==0 && this->M_dataFluid->isSemiImplicit())|| !this->M_dataFluid->isSemiImplicit())
            {

                vector_type lambdaFluid(this->M_interfaceMap, Unique);

                monolithicToInterface(lambdaFluid, disp);

                lambdaFluid *= M_dataFluid->timestep();

                this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface

                //                M_epetraWorldComm->Barrier();
                chronoFluid.start();

                M_meshMotion->iterate();

                //                this->transferMeshMotionOnFluid(M_meshMotion->disp(),
                //                                        this->veloFluidMesh());

                //                this->veloFluidMesh()    -= dispFluidMeshOld();
                //                this->veloFluidMesh()    *= 1./(M_dataFluid->timestep());

                // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
                // the value of the displacement for some points
                                *this->M_beta *= 0.;
                                *this->M_rhs *= 0.;

                vector_type const meshDispDiff( M_meshMotion->dispDiff(), Repeated );

                this->moveMesh(meshDispDiff);

                this->interpolateVelocity(meshDispDiff, *this->M_beta);

                *this->M_beta *= -1./M_dataFluid->timestep();
                vector_ptrtype fluid(new vector_type(this->M_uFESpace->map()));
                fluid->subset(*M_un, 0);
                *this->M_beta += *fluid/*M_un*/;

                double alpha = 1./M_dataFluid->timestep();

                M_fluid->updateSystem(alpha,* this->M_beta, *this->M_rhs );//here it assembles the fluid matrices

                // coupling matrices assebling

                *this->M_rhs               += M_fluid->matrMass()*M_bdf->time_der( M_dataFluid->timestep() );
                matrix_ptrtype monolithicMatrix;
                monolithicMatrix.reset(new matrix_type(M_monolithicMap/*, this->M_fluid->getMeanNumEntries()*/));

                //*monolithicMatrix += *this->M_fluid->getBigMatrixPtr();
                M_fluid->updateStab( *monolithicMatrix);//applies the stabilization terms
                M_fluid->getFluidMatrix( *monolithicMatrix);

                coupling(monolithicMatrix, this->M_rhs);

                /*
                if(firstIter)
                    {
                        M_solid->buildSystem( *monolithicMatrix);
                        firstIter=false;
                    }
                else
                    M_solid->updateMatrix( *monolithicMatrix);
                */

                if(firstIter)
                    {
                        M_solid->buildSystem();
                        M_solid->rescaleMatrices();
                        firstIter=false;
                    }

                this->M_solid->updateMatrix( *monolithicMatrix);


                this->setDispSolid();

                this->M_solid->updateSystem(*this->M_rhs);

                this->M_solid->evalResidual( *M_BCh_u, *M_BCh_d, disp, res/*,  *this->M_fluid->getBigMatrixPtr() */);
            }

        vector_ptrtype residual(new vector_type(this->M_monolithicMap, Unique));
        this->M_solid->evalResidual( disp, *residual );
        res = *residual;

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
            bigPrecPtr.reset(new matrix_type(M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
            this->M_fluid->setBlockPreconditioner(bigPrecPtr);
                //            *bigPrecPtr += *this->M_fluid->getBlockPrecPtr();
            this->M_fluid->updateStab(*bigPrecPtr);
            Real entry(1.0);
            this->M_solid->setBlockPreconditioner(bigPrecPtr);
            addDiagonalEntries(entry,bigPrecPtr);
            if(!M_isDiagonalBlockPrec)
                coupling(bigPrecPtr);
            bigPrecPtr->GlobalAssemble();
        }
    else
        {
            if(M_isDiagonalBlockPrec)
                {
                    bigPrecPtr.reset(new matrix_type(M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
                    this->M_fluid->updateStab( *bigPrecPtr);//applies the stabilization terms
                    this->M_fluid->getFluidMatrix( *bigPrecPtr);
                    this->M_solid->setBlockPreconditioner(bigPrecPtr);
                    Real entry(1.0);
                    addDiagonalEntries(entry,bigPrecPtr);
                    bigPrecPtr->GlobalAssemble();
                }
            else
                {
                    bigPrecPtr.reset(new matrix_type(M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
                    this->M_solid->setFullPreconditioner(bigPrecPtr);
                    Real entry(1.0);
                    addDiagonalEntries(entry,bigPrecPtr);
                    bigPrecPtr->GlobalAssemble();
                }
        }
    this->M_solid->iterateMonolithic(*(const_cast<vector_type*>(&_res))/*just to avoid the const type*/, _step, bigPrecPtr);

    setDispSolid();// to be done just if not semiImplicit (done once yet in updateSystem)
    this->M_solid->updateStuff();
    std::cout << "done." << std::endl;

}

void
Monolithic::iterateMesh(const vector_type& disp)
{
                vector_type lambdaFluid(this->M_interfaceMap, Unique);

                monolithicToInterface(lambdaFluid, disp);

                lambdaFluid *= M_dataFluid->timestep();

                this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface

                //                M_epetraWorldComm->Barrier();

                M_meshMotion->iterate();

                //                this->transferMeshMotionOnFluid(M_meshMotion->disp(),
                //                                        this->veloFluidMesh());

                //                this->veloFluidMesh()    -= dispFluidMeshOld();
                //                this->veloFluidMesh()    *= 1./(M_dataFluid->timestep());

                // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
                // the value of the displacement for some points
                //                vector_type const meshDispDiff( M_meshMotion->dispDiff(), Repeated );

                //                this->moveMesh(meshDispDiff);
}



namespace
{
FSIOperator* createM(){ return new Monolithic(); }
}
static bool reg = FSIFactory::instance().registerProduct( "monolithic", &createM );

}
