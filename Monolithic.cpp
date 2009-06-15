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

#include <lifemc/lifesolver/Monolithic.hpp>
#include "EpetraExt_MatrixMatrix.h"
//#include <life/lifesolver/reducedLinFluid.hpp>

namespace LifeV
{
// Constructors
Monolithic::Monolithic():
    super(),
    M_monolithicMap(),
    M_couplingMatrix(),
    M_interface(0),
    M_interfaceMap(),
    M_beta(),
    M_monolithicMatrix(),
    M_bigPrecPtr(),
    M_DDBlockPrec(),
    M_rhsFull(),
    //end of protected args
    M_updateEvery(0),
    M_numerationInterface(),
    M_rhsShapeDerivatives(),
    M_rhsNew(),
    M_fullMonolithic(false),
    M_entry(0.),
    M_robinCoupling(),
    //    M_robinCouplingInv(),
    M_alphaf(0.),
    M_alphas(0.),
    //    M_isDiagonalBlockPrec(false),
    M_diagonalScale(false),
    M_offset(0),
    M_solidAndFluidDim(0),
    M_reusePrec(true),
    M_resetPrec(true),
    M_maxIterSolver(-1),
    M_linearSolver()
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

            M_monolithicMap.reset(new EpetraMap(M_uFESpace->map()));
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
            couplingVector.reserve((int)(M_interfaceMap.getMap(Unique)->NumMyElements()));

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
            M_beta.reset  (new vector_type(/*M_monolithicMap*/M_uFESpace->map()));

            M_offset = M_uFESpace->dof().numTotalDof()*nDimensions +  M_pFESpace->dof().numTotalDof();
            M_solidAndFluidDim= M_offset + M_dFESpace->dof().numTotalDof()*nDimensions;
            M_BCh_d->setOffset(M_offset);
            if(!M_fullMonolithic)
                {
                    M_meshMotion.reset(new FSIOperator::meshmotion_raw_type(*M_mmFESpace,
                                                                            *M_epetraComm));
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
            M_rhsFull.reset(new vector_type(*this->M_monolithicMap));
            M_un.reset (new vector_type(*this->M_monolithicMap));
            M_beta.reset  (new vector_type(/*M_monolithicMap*/M_uFESpace->map()));

            M_solid.reset(new FSIOperator::solid_raw_type(dataSolid(),
                                                         *M_dFESpace,
                                                         *M_epetraComm,
                                                         *M_monolithicMap,
                                                         M_offset
                                                         ));

//             if (isLinearSolid())// to be implemented with the offset
//                 M_solidLin.reset(new FSIOperator::solidlin_raw_type(dataSolid(),
//                                                                    *M_dFESpace,
//
//                                                      *M_epetraComm));
                }
}//end setup

void
Monolithic::setDataFromGetPot( GetPot const& data_file )
{
    super::setDataFromGetPot(data_file);

    M_updateEvery = data_file("problem/updateEvery", 0);
    this->M_dataFluid->setSemiImplicit( data_file("problem/semiImplicit", true) );
    this->M_dataFluid->setUseShapeDerivatives( data_file("fluid/useShapeDerivatives", false) );
    M_DDBlockPrec = data_file( "interface/DDBlockPrec",  0 );
    //M_splitPrec = data_file( "interface/splitPrec",  false );
    M_fullMonolithic  = !(M_method.compare("fullMonolithic"));
    M_diagonalScale           = data_file( "solid/prec/diagonalScaling",  false );
    M_entry           = data_file( "solid/prec/entry",  0. );
    M_alphaf           = data_file( "interface/alphaf",  0.5 );
    M_alphas           = data_file( "interface/alphas",  0.5 );
}

void
Monolithic::buildSystem()
{
    M_couplingMatrix.reset(new matrix_type(*M_monolithicMap/*, this->M_fluid->getMeanNumEntries()*/));
    this->couplingMatrix(M_couplingMatrix);
    M_couplingMatrix->GlobalAssemble();

    M_solid->buildSystem();
    M_solid->rescaleMatrices();
    updateCoupling(*M_couplingMatrix);
}

void
Monolithic::couplingMatrix(matrix_ptrtype & bigMatrix, bool solidCoupling) // not working with non-matching grids
{

    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;

    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
                    if(solidCoupling)
                        {
                            bigMatrix->set_mat_inc( M_offset + ITrow->second-1 + dim* M_dFESpace->dof().numTotalDof(),(int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (1.) );//right low
                        }

                            bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (M_offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), (-1.*M_solid->rescaleFactor()/*/M_dataFluid->timestep()*/));//low right
                    bigMatrix->set_mat_inc( ITrow->first-1 + dim* M_uFESpace->dof().numTotalDof(), (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (-1.) );//right up
                    bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), 1.);//low left
                    //                    bigMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim , (int)(*this->M_numerationInterface)[ITrow->second /*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, 0.0);
                }
        }
    }
}


void
Monolithic::robinCoupling(matrix_ptrtype & IdentityMatrix, Real& alphaf, Real& alphas/*, bool inverse*/) // not working with non-matching grids
{
    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    //Real* entry(new Real(1.));

    addDiagonalEntries(1., IdentityMatrix, *M_monolithicMap);
    //addDiagonalEntries(alphaf, IdentityMatrix, *M_solidInterfaceMap, M_offset, true);
    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
// //                     if(inverse)
// //                         {//pay attention to the indices!! alphas is right low
// //                         IdentityMatrix->getEpetraMatrix().ReplaceGlobalValues( M_offset + ITrow->second + dim* M_dFESpace->dof().numTotalDof(),M_offset + ITrow->second + dim* M_dFESpace->dof().numTotalDof(), (-alphaf) );//diagonal up
// //                         IdentityMatrix->getEpetraMatrix().ReplaceGlobalValues( ITrow->first + dim* M_uFESpace->dof().numTotalDof(),ITrow->first + dim* M_uFESpace->dof().numTotalDof() , (1-alphas) );//diagonal low
// //                         }
                    //IdentityMatrix->set_mat_inc( ITrow->first-1 + dim* M_uFESpace->dof().numTotalDof(),M_offset + ITrow->second-1 + dim* M_dFESpace->dof().numTotalDof() , (alphaf));//right up
//                     IdentityMatrix->set_mat_inc( ITrow->second-1 + dim* M_dFESpace->dof().numTotalDof()+M_offset, (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), alphas);//left low
//                     //                    IdentityMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim , (int)(*this->M_numerationInterface)[ITrow->second /*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, 0.0);
                    //IdentityMatrix->set_mat_inc( M_offset + ITrow->second-1 + dim* M_dFESpace->dof().numTotalDof(),(int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (1.) );//right low
                    IdentityMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (M_offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), alphas/*(1.)*/);//low right
                    IdentityMatrix->set_mat_inc( ITrow->first-1 + dim* M_uFESpace->dof().numTotalDof(), (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, alphaf/*(1.)*/ );//right up
                    IdentityMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), alphas);//low left
                    //int* index(new int((int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] + dim*M_interface + M_solidAndFluidDim));
                    //IdentityMatrix->getEpetraMatrix().ReplaceGlobalValues((*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] + dim*M_interface + M_solidAndFluidDim,1, entry,  index);//(5, 5)
                }
        }
    }
}
void
Monolithic::robinCouplingInv(matrix_ptrtype & IdentityMatrix, Real& alphaf, Real& alphas/*, bool inverse*/) // not working with non-matching grids
{
    Real factor = 1/(alphas*alphaf+1);
    Real y=alphaf*factor;
    Real l=alphas*factor;
    Real s=alphas*alphaf*factor;
    Real* entry(new Real(-y));

    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;
    addDiagonalEntries(1., IdentityMatrix, *M_monolithicMap);
    for(short i=0; i<nDimensions; ++i)
        addDiagonalEntries(factor, IdentityMatrix, locDofMap, i*M_uFESpace->dof().numTotalDof(), true);//(2,2)
    addDiagonalEntries(l, IdentityMatrix, *M_solidInterfaceMap, M_offset, true);//(4,4)

    for(UInt dim = 0; dim < nDimensions; ++dim)
    {
        for( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                {
                    IdentityMatrix->set_mat_inc( M_offset + ITrow->second-1 + dim* M_dFESpace->dof().numTotalDof(),(int)(*this->M_numerationInterface)[ITrow->second ] - 1 + dim*M_interface + M_solidAndFluidDim, factor );//(4,5)
                    IdentityMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (M_offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), factor);//(5,4)
                    IdentityMatrix->set_mat_inc( ITrow->first-1 + dim* M_uFESpace->dof().numTotalDof(), (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, y );//(2,5)
                    IdentityMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), s);//(5,2)
                    //IdentityMatrix->set_mat_inc( (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim , (int)(*this->M_numerationInterface)[ITrow->second /*+ dim*solidDim*/ ] - 1 + dim*M_interface + M_solidAndFluidDim, -y); //(5,5)
                    IdentityMatrix->set_mat_inc( (M_offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), -l);//(4,2)
                    IdentityMatrix->set_mat_inc( (ITrow->first)-1 + dim* M_uFESpace->dof().numTotalDof(), (M_offset + ITrow->second)-1 + dim* M_dFESpace->dof().numTotalDof(), -factor);//(2,4)
                    int* index(new int((int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] + dim*M_interface + M_solidAndFluidDim));
                    IdentityMatrix->getEpetraMatrix().ReplaceGlobalValues((*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] + dim*M_interface + M_solidAndFluidDim,1, entry,  index);//(5,5)
                }
        }
    }
}

void Monolithic::zeroBlock( matrix_ptrtype matrixPtr,vector_type& colNumeration , const std::map<ID, ID>& map, UInt rowOffset, UInt colOffset)
{//to improve
    std::map<ID, ID> const& Map = M_dofStructureToHarmonicExtension->locDofMap();
    Real* entry(new Real(0.));
    std::map<ID, ID>::const_iterator ITrow;
    int err;
    if(colOffset>=M_offset && colOffset<M_solidAndFluidDim && rowOffset<M_solidAndFluidDim)
        {
            for( ITrow=map.begin(); ITrow != map.end(); ++ITrow)// scalable loops,
                {
                    //                    if(colNumeration.getMap().getMap(Unique)->LID(ITrow->second)>=0)
                     if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                        {
                            int index=(colOffset + ITrow->second);// second is column
                            err=matrixPtr->getEpetraMatrix().ReplaceGlobalValues(rowOffset + ITrow->first,1, entry ,&index);
                        }
                }
            return;
        }
    if(rowOffset>=M_offset && rowOffset<M_solidAndFluidDim && colOffset<M_solidAndFluidDim)
        {
            for( ITrow=map.begin(); ITrow != map.end(); ++ITrow)// scalable loops,
                {
                //                    if(colNumeration.getMap().getMap(Unique)->LID(ITrow->second)>=0)
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                        {
                            int index2= (colOffset + ITrow->first);
                            err=matrixPtr->getEpetraMatrix().ReplaceGlobalValues(rowOffset + ITrow->second,1, entry ,&index2);
                        }
                }
            return;
        }

    if(colOffset>=M_solidAndFluidDim && rowOffset < M_offset)
        {
            for( ITrow=map.begin(); ITrow != map.end(); ++ITrow)// scalable loops,
                {
                    //if(colNumeration.getMap().getMap(Unique)->LID(ITrow->second)>=0)
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                        {
                            int* index(new int(colOffset + colNumeration[ITrow->second]));
                           err= matrixPtr->getEpetraMatrix().ReplaceGlobalValues(rowOffset + ITrow->first,1, entry ,index);

                        }
                }
            return;
        }

    if((colOffset>=M_solidAndFluidDim && rowOffset >= M_offset)/*||(colOffset>=M_offset && rowOffset >= M_solidAndFluidDim)*/)
        {
            for( ITrow=map.begin(); ITrow != map.end(); ++ITrow)// scalable loops,
                {
                    //if(colNumeration.getMap().getMap(Unique)->LID(ITrow->second)>=0)
            if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                        {
                            int* index2(new int(colOffset + colNumeration[ITrow->second]));
                            err=matrixPtr->getEpetraMatrix().ReplaceGlobalValues( rowOffset + ITrow->second,1, entry ,index2);
                        }
            if (err != 0)
                {
                    M_solid->getDisplayer().leaderPrint("error ", err);
                    M_solid->getDisplayer().leaderPrint("Go to fill your matrix!!");
                }
                }
            return;
        }

}

void
Monolithic::couplingRhs(vector_ptrtype rhs, vector_ptrtype un) // not working with non-matching grids
{
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
                        (*rhs)[  (int)(*this->M_numerationInterface)[ITrow->second/*+ dim*solidDim*/ ] + dim*M_interface +M_solidAndFluidDim ] = -lambda( ITrow->second + dim*M_dFESpace->dof().numTotalDof()/*M_dFESpace->dof().numTotalDof()*/ )*(M_solid->rescaleFactor())/*/M_dataFluid->timestep()*/;
                    //                    std::cout<<(int)(*this->M_numerationInterface)[ITrow->second/*+dim*solidDim*/] - 1 + dim*M_interface<<std::endl;

                }
        }
    }
}


void
Monolithic::addDiagonalEntries(Real entry, matrix_ptrtype bigMatrix, EpetraMap& Map, UInt offset, bool replace)
{
    if(!replace)
        for(UInt i=0 ; i<Map.getMap(Repeated)->NumMyElements(); ++i)//num from 1
            {
                bigMatrix->set_mat_inc(  offset + Map.getMap(Repeated)->GID(i)-1 ,   offset + Map.getMap(Repeated)->GID(i)-1, entry);
            }
    else
        {
            for(UInt i=0 ; i<Map.getMap(Unique)->NumMyElements(); ++i)
                {
                    //diagonal[Map.getMap(Repeated)->GID(i)+offset]=entry;
                    int* index(new int(offset + Map.getMap(Unique)->GID(i)));
                    bigMatrix->getEpetraMatrix().ReplaceGlobalValues(offset + Map.getMap(Repeated)->GID(i), 1, &entry, index );
                }
        }
}
void
Monolithic::addDiagonalEntries(Real entry,  matrix_ptrtype bigMatrix, std::map<ID, ID> const& Map, UInt offset, bool replace)
{
    std::map<ID, ID>::const_iterator ITrow;
    if(!replace)
        for(ITrow=Map.begin() ; ITrow != Map.end(); ++ITrow)
            {
                if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                bigMatrix->set_mat_inc(  offset + (*M_numerationInterface)[ITrow->second]-1 ,   offset + (*M_numerationInterface)[ITrow->second]-1, entry);
            }
    else
        {
            for(ITrow=Map.begin() ; ITrow != Map.end(); ++ITrow)
                {
                    if(M_interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/) >= 0 )//to avoid repeated stuff
                    {
                        int* index(new int(offset + (*M_numerationInterface)[ITrow->second]));
                        int err = bigMatrix->getEpetraMatrix().ReplaceGlobalValues(offset + (*M_numerationInterface)[ITrow->second], 1, &entry, index );
                        if (err != 0)
                            {
                                M_solid->getDisplayer().leaderPrint("error ", err);
                                M_solid->getDisplayer().leaderPrint("Go to fill your matrix!!");
                            }
                    }
                }
        }
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
    *M_rhsFull*=0;
    this->M_fluid->resetStab();
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
    //if (dispFluid.getMaptype() == Unique)
           /*        {
            vector_type  dispFluidRep(dispFluid.getMap(), Repeated);
            monolithicToFluid(disp, dispFluidRep);
            dispFluid = dispFluidRep;
            return;
            }*/

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

       EpetraMap subMap(*disp.getMap().getMap(Unique), M_offset,disp.getMap().getMap(Unique)->NumGlobalElements() );
       vector_type subDisp(subMap, Unique);
       subDisp.subset(disp, M_offset);
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
            dispSolid = dispSolidUn;
            return;
        }

    ASSERT( &disp != &dispSolid, "Calling monolithicToSolid with the same two arguments");

    dispSolid *= 0.;
    UInt k = 0;
    UInt N = M_dFESpace->dof().numTotalDof();

    for( UInt dim = 0; dim<nDimensions; ++dim)
        {
            for( k=1 ; k <= N ; ++k)// since map iterators start from one, operator [] subtracts one!
                {
                    if(M_monolithicMap->getMap(Unique)->LID(k+N*dim+M_offset) >= 0 )
                        {
                            dispSolid[k + N*dim + M_offset] = disp(k + N*dim + M_offset);
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
                          const int          iter )
{
            //            eval(disp, iter);

    Chrono chronoFluid, chronoSolid, chronoInterface;


    if((iter==0)|| !this->M_dataFluid->isSemiImplicit())
        {
            setDispSolid(disp);
            //            vector_ptrtype meshDeltaDisp;//is the meshDispOld used to calculate deltaDisp
            vector_type lambdaFluid(this->M_interfaceMap, Unique);

            monolithicToInterface(lambdaFluid, disp);

            lambdaFluid *= (M_dataFluid->getTimeStep()*(M_solid->rescaleFactor()));//because of the matrix scaling

            this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface

            //                M_epetraWorldComm->Barrier();
            chronoFluid.start();
            M_meshMotion->iterate();
            M_meshMotion->updateDispDiff();
            //                this->transferMeshMotionOnFluid(M_meshMotion->disp(),
            //                                        this->veloFluidMesh());

            //                this->veloFluidMesh()    -= dispFluidMeshOld();
            //                this->veloFluidMesh()    *= 1./(M_dataFluid->timestep());

            // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
            // the value of the displacement for some points
            //if(iter==0)
            //{
                    //                    M_rhs.reset(new vector_type(*M_monolithicMap));
            //}
            M_beta.reset(new vector_type(M_uFESpace->map()));

            vector_type const meshDispDiff( M_meshMotion->dispDiff(), Repeated );

            this->moveMesh(meshDispDiff);

            this->interpolateVelocity(meshDispDiff, *this->M_beta);

            double alpha = 1./M_dataFluid->getTimeStep();

            *this->M_beta *= -alpha;
            vector_ptrtype fluid(new vector_type(this->M_uFESpace->map()));
            fluid->subset(*M_un, 0);
            *this->M_beta += *fluid/*M_un*/;

            //            M_beta->GlobalAssemble();

            M_fluid->updateSystem(alpha,*this->M_beta, *this->M_rhs/*, M_monolithicMatrix*/ );//here it assembles the fluid matrices
            M_monolithicMatrix.reset(new matrix_type(*M_monolithicMap));

            // coupling matrices assebling

            M_fluid->getFluidMatrix( *M_monolithicMatrix );

            M_fluid->updateStab( *M_monolithicMatrix);//applies the stabilization terms

            this->updateMatrix( *M_monolithicMatrix);
            M_monolithicMatrix->GlobalAssemble();
            //M_epetraWorldComm->Barrier();
            if(iter==0)
                {
                    M_nbEval = 0; // new time step
                    M_resetPrec=true;
                    *this->M_rhs               += M_fluid->matrMass()*M_bdf->time_der( M_dataFluid->getTimeStep() );
                    couplingRhs(this->M_rhs, M_un);
                    this->M_solid->updateStuff();
                    updateSolidSystem(this->M_rhs);
                    //M_rhs.reset(new vector_type(*M_monolithicMap));
                }
            M_rhsFull.reset(new vector_type(*M_rhs));
            evalResidual( *M_BCh_u, *M_BCh_d, disp, M_rhsFull, res, M_diagonalScale);

            if(M_DDBlockPrec>=2 && M_DDBlockPrec<5)
            {
                if(!M_robinCoupling.get())
                    {
                        M_robinCoupling.reset( new matrix_type(*M_monolithicMap));
                        robinCoupling(M_robinCoupling, M_alphaf, M_alphas);
                        M_robinCoupling->GlobalAssemble();
                    }
                this->applyPreconditioner(M_robinCoupling);
                //this->M_solid->evalResidual( disp, res, false);
            }

            if(M_DDBlockPrec<5)
                M_bigPrecPtr.reset(new matrix_type(*M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));

            M_nbEval++ ;
        }
    //else
    //{
    evalResidual( disp,  M_rhsFull, res, false);
            //}
}

void Monolithic::shapeDerivatives(vector_ptrtype rhs, vector_ptrtype meshDeltaDisp, const vector_type& sol)
{
    double alpha = 1./M_dataFluid->getTimeStep();
    vector_type derVeloFluidMsh(*meshDeltaDisp, Repeated);
    derVeloFluidMsh *= alpha;
    vector_ptrtype rhsNew(new vector_type(*M_monolithicMap));
    vector_type un(M_uFESpace->map()/*+M_pFESpace->map()*/);
    vector_type uk(M_uFESpace->map()+M_pFESpace->map());
    vector_type meshVel(M_meshMotion->dispDiff(), Repeated);
    un.subset(*this->M_un, 0);
    uk.subset(sol, 0);
    //    if(!M_rhsNew.get())
    //    M_rhsNew.reset(new vector_type(M_uFESpace->map()+M_pFESpace->map()));//useless
    vector_type dvfm(M_uFESpace->map()/*+M_pFESpace->map()*/, Repeated);
    vector_type vfm(M_uFESpace->map()/*+M_pFESpace->map()*/, Repeated);
    vector_type vmdisp(M_uFESpace->map()/*+M_pFESpace->map()*/, Repeated);
    this->transferMeshMotionOnFluid(derVeloFluidMsh, dvfm);
    this->transferMeshMotionOnFluid(*meshDeltaDisp, vmdisp);
    this->transferMeshMotionOnFluid(meshVel, vfm);

    //    rhs->spy("rhsSd");
    //    meshDeltaDisp->spy("deltaDispSd");
    //    sol.spy("solSd");
    //    un.spy("unSd");
    //    uk.spy("ukSd");
    //    vmdisp.spy("deltax");
    //    vfm.spy("velFluidMesh");
    //    dvfm.spy("derVelFluidMesh");

    M_fluid->updateLinearSystem(*M_monolithicMatrix,
                                0.,
                                un,//un
                                uk,//uk
                                vmdisp,//unknown x.
                                //Repeated
                                vfm, //(xk-xn)/dt //Repeated
                                dvfm, // (x)/dt //Repeated
                                *rhsNew// useless
                                );
    //    M_rhsShapeDerivatives.reset(new vector_type(*this->M_rhs)) ;
    *rhs = M_fluid->rhsLinNoBC();// Import
    //    rhs->spy("rhsShape");
    std::cout<<"rhs lin no bc : "<<M_fluid->rhsLinNoBC().NormInf()<<std::endl;
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
evalResidual( fluid_bchandler_raw_type& bchFluid, solid_bchandler_raw_type& bchSolid, const vector_type& sol, vector_ptrtype& rhs, vector_type& res, bool diagonalScaling, matrix_ptrtype preconditioner)
{
    //M_solid->getMatrixPtr()->GlobalAssemble();
    M_monolithicMatrix->GlobalAssemble();
    matrix_ptrtype tmpMatPtr(new matrix_type(*M_monolithicMatrix));
    tmpMatPtr->GlobalAssemble();
    //preconditioner->GlobalAssemble();
    //M_solid->getMatrixPtr().reset(new matrix_type(M_solid->getMap()));
    M_monolithicMatrix.reset(new matrix_type(M_solid->getMap()));
    int err = EpetraExt::MatrixMatrix::Multiply( preconditioner->getEpetraMatrix(), false, tmpMatPtr->getEpetraMatrix(), false, M_monolithicMatrix->getEpetraMatrix());
    *rhs = (*preconditioner)*(*rhs);
    evalResidual( bchFluid, bchSolid, sol, rhs, res, diagonalScaling);
}

void Monolithic::
evalResidual( fluid_bchandler_raw_type& bchFluid, solid_bchandler_raw_type& bchSolid, const vector_type& sol, vector_ptrtype& rhs, vector_type& res, bool diagonalScaling)
{
    vector_type rhsFullSolid(*rhs, Unique); // ignoring non-local entries, Otherwise they are summed up lately
    //    applyBoundaryConditions( *M_solid->getMatrixPtr(), rhsFullSolid, bchSolid, M_solid->offset() );
        if(M_solid->offset())
            bchSolid.setOffset(M_solid->offset());
        if ( !bchSolid.bdUpdateDone() )
            bchSolid.bdUpdate( *M_dFESpace->mesh(), M_dFESpace->feBd(), M_dFESpace->dof() );

        bcManage( *M_monolithicMatrix, rhsFullSolid, *M_dFESpace->mesh(), M_dFESpace->dof(), bchSolid, M_dFESpace->feBd(), 1.,
                  dataSolid().getTime() );

    // matrix is GlobalAssembled by  bcManage

    if ( !bchFluid.bdUpdateDone() )
        bchFluid.bdUpdate( *M_uFESpace->mesh(), M_uFESpace->feBd(), M_uFESpace->dof() );

    //vector_type rhsFull(rhsFullSolid, Unique);
    vector_type rhsFull(rhsFullSolid);
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

    M_solid->getDisplayer().leaderPrint("solveJac: NormInf res ", _res.NormInf());
    //_step *= 0.;

     M_solid->getDisplayer().leaderPrint("Solving Jacobian system... \n" );
    //    M_epetraWorldComm->Barrier();

     //    matrix_ptrtype bigPrecPtr;

    if(this->M_fluid->getIsDiagonalBlockPrec()==true /*M_diagonalBlockPreconditioner*/)
        {
            M_bigPrecPtr.reset(new matrix_type(*M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
            this->M_fluid->setBlockPreconditioner(M_bigPrecPtr);
                //            *bigPrecPtr += *this->M_fluid->getBlockPrecPtr();
            this->M_fluid->updateStab(*M_bigPrecPtr);
            Real entry(1.0);
            this->setBlockPreconditioner(M_bigPrecPtr);
            addDiagonalEntries(entry,M_bigPrecPtr, M_interfaceMap);
            //if(!M_isDiagonalBlockPrec)
            couplingMatrix(M_bigPrecPtr);
            M_bigPrecPtr->GlobalAssemble();
        }
    else
        {
            if( !M_reusePrec || M_resetPrec )
                        switch(M_DDBlockPrec)
                            {
                            case 1:
                                {
                                    //M_bigPrecPtr.reset(new matrix_type(*M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
                                    this->M_fluid->getFluidMatrix( *M_bigPrecPtr);
                                    this->M_fluid->updateStab( *M_bigPrecPtr);//applies the stabilization terms
                                    if(!M_bigPrecPtr->getEpetraMatrix().Filled())
                                    {
                                        couplingMatrix(M_bigPrecPtr, false);
                                        addDiagonalEntries(M_entry,M_bigPrecPtr, M_interfaceMap, M_solidAndFluidDim);
                                    }
                                    else
                                        this->M_solid->getDisplayer().leaderPrint("ERROR: probably the tolerance fixed for Newton is too low. \n");
                                    this->setBlockPreconditioner(M_bigPrecPtr);
                                    //this->M_solid->setFullPreconditioner(M_bigPrecPtr);
                                    //std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
                                    //for(short i=0; i<nDimensions; ++i)
                                    //    zeroBlock(M_bigPrecPtr, *M_numerationInterface, locDofMap, M_offset + i*M_dFESpace->dof().numTotalDof(), i*M_interface+M_solidAndFluidDim);//(4, 5)
                                    ////Real entry(0.0);

                                    M_bigPrecPtr->GlobalAssemble();}
                                break;
                            case 2:
                                {
                                    //                                 if(false && !M_robinCoupling.get())
                                    //                                 {
                                    //                                     ////Real factor=1/(1-M_alphaf*M_alphas);
                                    //                                     M_robinCoupling.reset( new matrix_type(*M_monolithicMap));
                                    //                                     robinCoupling(M_robinCoupling, M_alphaf, M_alphas);
                                    //                                     ////addDiagonalEntries((Real)1.0, M_robinCoupling, *M_monolithicMap, (UInt)0);
                                    //                                     M_robinCoupling->GlobalAssemble();
                                    //                                     M_robinCouplingInv.reset( new matrix_type(*M_monolithicMap/**M_robinCoupling*/));
                                    //                                     robinCouplingInv(M_robinCouplingInv, M_alphaf, M_alphas);

                                    //                                     ////addDiagonalEntries((Real)1.0/factor, M_robinCouplingInv, *M_fluidInterfaceMap, (UInt)0, true);
                                    //                                     ////addDiagonalEntries((Real)1.0/factor, M_robinCouplingInv, *M_solidInterfaceMap, M_offset, true);
                                    //                                     M_robinCouplingInv->GlobalAssemble();
                                    //                                     //TEST:
                                    //                                     matrix_ptrtype identity(new matrix_type(*M_monolithicMap));
                                    //                                     //addDiagonalEntries(1., identity, *M_monolithicMap, 0.);
                                    //                                     //identity->GlobalAssemble();
                                    //                                     //identity->spy("id");
                                    //                                     EpetraExt::MatrixMatrix::Multiply( M_robinCouplingInv->getEpetraMatrix(), false, M_robinCoupling->getEpetraMatrix(), false, identity->getEpetraMatrix());
                                    //                                     identity->spy("identity");
                                    //                                     //END TEST
                                    //                                 }
                                    //                             M_robinCoupling->spy("robinCoupling");
                                    //                             M_robinCouplingInv->spy("robinCouplingInv");
                                    //                             matrix_ptrtype tmpPrecPtr(new matrix_type(*M_monolithicMap));
                                    //                             M_bigPrecPtr.reset(new matrix_type(*M_monolithicMap));
                                    //                             this->M_solid->setFullPreconditioner(tmpPrecPtr);
                                    //                             tmpPrecPtr->getEpetraMatrix().FillComplete();
                                    //                             int err = EpetraExt::MatrixMatrix::Multiply( M_robinCoupling->getEpetraMatrix(), false, tmpPrecPtr->getEpetraMatrix(), false, M_bigPrecPtr->getEpetraMatrix()/*, false*/);
                                    //                             M_bigPrecPtr->spy("beforZeroBlock");

                                    //                             tmpPrecPtr.reset(new matrix_type(*M_bigPrecPtr));
                                    //                             //tmpPrecPtr->zeroBlock( 0, M_offset, M_offset, M_solidAndFluidDim);
                                    //                             std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
                                    //                             for(short i=0; i<nDimensions; ++i)
                                    //                                 {
                                    //                                     //zeroBlock(tmpPrecPtr, *M_numerationInterface, locDofMap, i*M_uFESpace->dof().numTotalDof(), i*M_dFESpace->dof().numTotalDof()+M_offset);
                                    //                                     //zeroBlock(tmpPrecPtr, *M_numerationInterface, locDofMap, i*M_dFESpace->dof().numTotalDof()+M_offset, i*M_uFESpace->dof().numTotalDof());
                                    //                                     ////M_bigPrecPtr->zeroBlock(locDofMap, i*M_uFESpace->dof().numTotalDof(), i*M_dFESpace->dof().numTotalDof()+M_offset);
                                    //                                     //zeroBlock(tmpPrecPtr, *M_numerationInterface, locDofMap, i*M_uFESpace->dof().numTotalDof(), i*M_interface+M_solidAndFluidDim);
                                    //                                     //zeroBlock(tmpPrecPtr, *M_numerationInterface, locDofMap, M_offset + i*M_dFESpace->dof().numTotalDof(), i*M_interface+M_solidAndFluidDim);
                                    //                                     ////M_bigPrecPtr->zeroBlock(*M_numerationInterface, locDofMap, i*M_uFESpace->dof().numTotalDof(), i*M_interface+M_solidAndFluidDim);
                                    //                                     //zeroBlock(tmpPrecPtr, *M_numerationInterface, locDofMap, M_solidAndFluidDim + i*M_interface, M_offset + i*M_dFESpace->dof().numTotalDof());

                                    //                                 }
                                    //                             tmpPrecPtr->GlobalAssemble();
                                    //                             M_bigPrecPtr->GlobalAssemble();
                                    //                             ////tmpPrecPtr->getEpetraMatrix().FillComplete();
                                    //                             tmpPrecPtr->spy("afterZeroBlock");
                                    //                             M_bigPrecPtr.reset(new matrix_type(*tmpPrecPtr/**M_monolithicMap*/));
                                    //                             ////M_bigPrecPtr.reset(new matrix_type(*tmpPrecPtr));
                                    //                             //err = EpetraExt::MatrixMatrix::Multiply( M_robinCouplingInv->getEpetraMatrix(), false, tmpPrecPtr->getEpetraMatrix(), false, M_bigPrecPtr->getEpetraMatrix());
                                    //                             M_bigPrecPtr->GlobalAssemble();
                                    //                             //M_solid->applyPreconditioner(*(const_cast<vector_type*>(&_res)), M_robinCoupling);

                                    //M_bigPrecPtr.reset(new matrix_type(*M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
                                    this->setFullPreconditioner(M_bigPrecPtr);
                                    //Real entry(0./*1.0e-5*/);
                                    //addDiagonalEntries(M_entry,M_bigPrecPtr, M_interfaceMap, M_solidAndFluidDim);
                                    M_bigPrecPtr->GlobalAssemble();
                                }
                                break;
                            case 3:
                                {if(false && !M_robinCoupling.get())
                                    {
                                        M_robinCoupling.reset( new matrix_type(*M_monolithicMap));
                                        robinCoupling(M_robinCoupling, M_alphaf, M_alphas);
                                        M_robinCoupling->GlobalAssemble();
                                    }
                                //matrix_ptrtype tmpPrecPtr(new matrix_type(*M_monolithicMap));
                                //M_bigPrecPtr.reset(new matrix_type(*M_monolithicMap));
                                this->setFullPreconditioner(M_bigPrecPtr/*tmpPrecPtr*/);


                                //this->M_fluid->setBlockPreconditioner(M_bigPrecPtr);
                                ////*M_bigPrecPtr += *this->M_fluid->getBlockPrecPtr();
                                //this->M_fluid->updateStab(*M_bigPrecPtr);
                                //Real entry(1.0);
                                //this->M_solid->setBlockPreconditioner(M_bigPrecPtr);
                                //addDiagonalEntries(M_entry, M_bigPrecPtr, M_interfaceMap, M_solidAndFluidDim);


                                //tmpPrecPtr->getEpetraMatrix().FillComplete();
                                //int err = EpetraExt::MatrixMatrix::Multiply( M_robinCoupling->getEpetraMatrix(), false, tmpPrecPtr->getEpetraMatrix(), false, M_bigPrecPtr->getEpetraMatrix()/*, false*/);
                                //M_bigPrecPtr->GlobalAssemble();
                                //M_bigPrecPtr->spy("beforZeroBlock");

                                //tmpPrecPtr.reset(new matrix_type(*M_bigPrecPtr));
                                //tmpPrecPtr->zeroBlock( 0, M_offset, M_offset, M_solidAndFluidDim);
                                std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
                                for(short i=0; i<nDimensions; ++i)
                                    {
                                        zeroBlock(M_bigPrecPtr/*tmpPrecPtr*/, *M_numerationInterface, locDofMap, i*M_uFESpace->dof().numTotalDof(), i*M_dFESpace->dof().numTotalDof()+M_offset); //(2, 4)
                                        zeroBlock(M_bigPrecPtr/*tmpPrecPtr*/, *M_numerationInterface, locDofMap, i*M_uFESpace->dof().numTotalDof(), i*M_interface+M_solidAndFluidDim);//(2, 5)
                                    }
                                //tmpPrecPtr->GlobalAssemble();
                                //M_bigPrecPtr->GlobalAssemble();
                                ////tmpPrecPtr->getEpetraMatrix().FillComplete();
                                //tmpPrecPtr->spy("afterZeroBlock");
                                //M_bigPrecPtr.reset(new matrix_type(*tmpPrecPtr/**M_monolithicMap*/));
                                //M_bigPrecPtr.reset(new matrix_type(*tmpPrecPtr));
                                //err = EpetraExt::MatrixMatrix::Multiply( M_robinCouplingInv->getEpetraMatrix(), false, tmpPrecPtr->getEpetraMatrix(), false, M_bigPrecPtr->getEpetraMatrix());
                                M_bigPrecPtr->GlobalAssemble();
                                //M_solid->applyPreconditioner(*(const_cast<vector_type*>(&_res)), M_robinCoupling);
                                }
                                break;
                            case 4:
                                {
                                    //M_bigPrecPtr.reset(new matrix_type(*M_monolithicMap));
                                    this->setFullPreconditioner(M_bigPrecPtr/*tmpPrecPtr*/);


                                    //                             this->M_fluid->updateStab(*M_bigPrecPtr);
                                    //                             this->M_fluid->getFluidMatrix( *M_bigPrecPtr);
                                    //                             this->M_solid->setBlockPreconditioner(M_bigPrecPtr);
                                    //                             addDiagonalEntries(M_entry, M_bigPrecPtr, M_interfaceMap, M_solidAndFluidDim);
                                    //                             couplingMatrix(M_bigPrecPtr, true);
                                    //                             M_bigPrecPtr->GlobalAssemble();
                                    //                             matrix_type tmpMatrix(*M_monolithicMap);
                                    //                             //tmpMatrix.GlobalAssemble();
                                    //                             EpetraExt::MatrixMatrix::Multiply( M_robinCoupling->getEpetraMatrix(), false, M_bigPrecPtr->getEpetraMatrix(), false, tmpMatrix.getEpetraMatrix());
                                    //                             M_bigPrecPtr.reset(new matrix_type(tmpMatrix));
                                    //                             M_bigPrecPtr->GlobalAssemble();
                                    M_epetraComm->Barrier();
                                    M_bigPrecPtr->GlobalAssemble();

                                    std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
                                    for(short i=0; i<nDimensions; ++i)
                                        {
                                            //zeroBlock(M_bigPrecPtr, *M_numerationInterface, locDofMap, i*M_uFESpace->dof().numTotalDof(), i*M_dFESpace->dof().numTotalDof()+M_offset);//(2, 4)
                                            //zeroBlock(M_bigPrecPtr, *M_numerationInterface, locDofMap, i*M_dFESpace->dof().numTotalDof()+M_offset, i*M_uFESpace->dof().numTotalDof());//(2, 5)
                                            //zeroBlock(M_bigPrecPtr, *M_numerationInterface, locDofMap, i*M_uFESpace->dof().numTotalDof(), i*M_interface+M_solidAndFluidDim);//(4, 2)
                                            zeroBlock(M_bigPrecPtr, *M_numerationInterface, locDofMap, M_offset + i*M_dFESpace->dof().numTotalDof(), i*M_interface+M_solidAndFluidDim);//(4, 5)
                                        }
                                }
                                break;
                            case 5:
                                {
                                }
                                break;
                            default:
                                {
                                    M_bigPrecPtr.reset(new matrix_type(*M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
                                    Real entry(0./*1.0e-5*/);
                                    for(short i=0; i<nDimensions; ++i)
                                        //addDiagonalEntries(M_entry,M_bigPrecPtr, M_dofStructureToHarmonicExtension->locDofMap(), M_solidAndFluidDim+i*M_interface, false);
                                    this->setFullPreconditioner(M_bigPrecPtr);
                                    //for(short i=0; i<nDimensions; ++i)
                                    //addDiagonalEntries(entry,M_bigPrecPtr, M_dofStructureToHarmonicExtension->locDofMap(), M_solidAndFluidDim+i*M_interface, true);
                                }
                                break;
                             }
        }
//             if(M_isDiagonalBlockPrec)
//                 {
//                     bigPrecPtr.reset(new matrix_type(M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
//                     this->M_fluid->updateStab( *bigPrecPtr);//applies the stabilization terms
//                     this->M_fluid->getFluidMatrix( *bigPrecPtr);
//                     this->setBlockPreconditioner(bigPrecPtr);
//                     Real entry(1.0);
//                     addDiagonalEntries(entry,bigPrecPtr);
//                     bigPrecPtr->GlobalAssemble();
//                 }
//             else
//                 {
//                     bigPrecPtr.reset(new matrix_type(M_monolithicMap/*, M_solid->getMatrixPtr()->getMeanNumEntries()*/));
//                     this->setFullPreconditioner(bigPrecPtr);
//                     Real entry(1.0);
//                     addDiagonalEntries(entry,bigPrecPtr);
//                    bigPrecPtr->GlobalAssemble();
//                }

    this->M_solid->getDisplayer().leaderPrint("preconditioner assembled \n");
    if(!M_dataFluid->useShapeDerivatives())
        setMatrix();
    this->iterateMonolithic(*(const_cast<vector_type*>(&_res))/*just to avoid the const type*/, _step, M_bigPrecPtr, M_linearSolver);

    //setDispSolid();// to be done just if not semiImplicit (done once yet in updateSystem)
    //this->M_solid->updateStuff();
    this->M_solid->getDisplayer().leaderPrint("done.\n");
}

void
Monolithic::iterateMesh(const vector_type& disp)
{
                vector_type lambdaFluid(this->M_interfaceMap, Unique);

                monolithicToInterface(lambdaFluid, disp);

                lambdaFluid *= (M_dataFluid->getTimeStep()*(M_solid->rescaleFactor()));

                this->setLambdaFluid(lambdaFluid); // it must be _disp restricted to the interface

                M_meshMotion->iterate();

}

void Monolithic::variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct)
{
    EpetraMap interfaceMap(*M_solidInterfaceMap);
    M_solidMeshPart.reset( new  partitionMesh< FSIOperator::mesh_type > (*M_dataSolid->mesh(), *M_epetraComm, M_solidInterfaceMap->getMap(Unique).get(), M_solidInterfaceMap->getMap(Repeated).get()));

    //M_solidMeshPart.reset( new  partitionMesh< mesh_type > (*M_dataSolid->mesh(), *M_epetraComm));

    M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(*M_solidMeshPart,
                                                       *refFE_struct,
                                                       *qR_struct,
                                                       *bdQr_struct,
                                                       3,
                                                       *M_epetraComm));
    // INITIALIZATION OF THE VARIABLES
    M_lambdaFluid.reset(new vector_type(*M_fluidInterfaceMap, Unique) );
    M_lambdaFluidRepeated.reset(new vector_type(*M_fluidInterfaceMap, Repeated) );

}

// void Monolithic::setUpSystem( GetPot const& data_file )
// {
//     super::setUpSystem(data_file);
//     M_solid->setAztecooPreconditioner(data_file, "solid/solver");
// }

void Monolithic::
updateCoupling(matrix_type couplingMatrix)//it works because it copies the matrix
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
    bigMatrixStokes += *M_couplingMatrix;
    //    bigMatrixStokes->GlobalAssemble();
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


void Monolithic::solidInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct)
{   // Monolitic: In the beginning I need a non-partitioned mesh. later we will do the partitioning
    M_dFESpace.reset(new FESpace<mesh_type, EpetraMap>(M_dataSolid->mesh(),
                                                       *refFE_struct,
                                                       *qR_struct,
                                                       *bdQr_struct,
                                                       3,
                                                       *M_epetraComm));
}

namespace
{
FSIOperator* createM(){ return new Monolithic(); }
}
static bool reg = FSIFactory::instance().registerProduct( "monolithic", &createM );

}
