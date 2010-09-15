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
    @date 08 Jun 2010

    A more detailed description of the file (if necessary)
 */

#include <BlockMatrix.hpp>

namespace LifeV {




int BlockMatrix::solveSystem( const vector_type& rhs, vector_type& step, solver_ptrtype& linearSolver)
{
    return linearSolver->solveSystem(rhs, step, M_globalMatrix);
}


void BlockMatrix::setDataFromGetPot(const GetPot& data, const std::string& section)
{
}


void BlockMatrix::coupler(map_shared_ptrtype& map,
                          const std::map<ID, ID>& locDofMap,
                          const vector_ptrtype& numerationInterface,
                          const Real& timeStep
                            )
{
    ASSERT(!M_coupling.get(), "coupler must not be called twice \n");
    M_coupling.reset(new matrix_type(*map));
    super::couplingMatrix( M_coupling,  M_couplingFlag, super::M_FESpace, super::M_offset, locDofMap, numerationInterface, timeStep);
}

void BlockMatrix::coupler(map_shared_ptrtype& map,
                          const std::map<ID, ID>& locDofMap,
                          const vector_ptrtype& numerationInterface,
                          const Real& timeStep,
                          UInt /*flag1*/
                            )
{
    super::couplingMatrix( M_coupling,  M_couplingFlag, super::M_FESpace, super::M_offset, locDofMap, numerationInterface, timeStep);
}


void BlockMatrix::push_back_matrix( const matrix_ptrtype& Mat, bool /*recompute*/)
{
    super::M_blocks.push_back(Mat);
}


void BlockMatrix::replace_matrix( const matrix_ptrtype& Mat, UInt index)
{
    super::M_blocks[index] = Mat;
}


void BlockMatrix::replace_precs( const epetra_operator_ptrtype& Mat, UInt index)
{
    assert(false);
}


void BlockMatrix::blockAssembling()
{
    M_coupling->GlobalAssemble();
    M_globalMatrix.reset(new matrix_type(M_coupling->getMap()));
    *M_globalMatrix += *M_coupling;
    for(UInt k=0; k<M_blocks.size(); ++k)
    {
        M_blocks[k]->GlobalAssemble();
        *M_globalMatrix += *M_blocks[k];
    }
}


void BlockMatrix::GlobalAssemble()
{
    M_globalMatrix->GlobalAssemble();
    //    M_globalMatrix->spy("Prec");
}


void BlockMatrix::applyPreconditioner(const matrix_ptrtype matrix, vector_ptrtype& rhsFull)
{
    this->applyPreconditioner(matrix, M_globalMatrix);
    *rhsFull = (*matrix)*(*rhsFull);
}


void BlockMatrix::applyPreconditioner( matrix_ptrtype robinCoupling, matrix_ptrtype prec, vector_ptrtype& rhs)
{
    applyPreconditioner(robinCoupling, prec);
    applyPreconditioner(robinCoupling, M_globalMatrix);
    (*rhs) = (*robinCoupling)*(*rhs);
}


void BlockMatrix::applyPreconditioner( const matrix_ptrtype prec, matrix_ptrtype& oper )
{
    matrix_type tmpMatrix(prec->getMap(), 1);
    EpetraExt::MatrixMatrix::Multiply( *prec->getMatrixPtr(),
                                       false,
                                       *oper->getMatrixPtr(),
                                       false,
                                       *tmpMatrix.getMatrixPtr());
    oper->swapCrsMatrix(tmpMatrix);
}


void BlockMatrix::createInterfaceMap( const EpetraMap& interfaceMap , const std::map<ID, ID>& locDofMap, const UInt subdomainMaxId,  const boost::shared_ptr<Epetra_Comm> epetraWorldComm )
{
    //std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;

    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    int* numInterfaceDof(new int[numtasks]);
    int pid=epetraWorldComm->MyPID();
    int numMyElements = interfaceMap.getMap(Unique)->NumMyElements();
    numInterfaceDof[pid]=numMyElements;
    EpetraMap subMap(*interfaceMap.getMap(Unique), (UInt)0, subdomainMaxId);

    M_numerationInterface.reset(new vector_type(subMap,Unique));
    //should be an int vector instead of double
    //                    M_numerationInterfaceInt.reset(new Epetra_IntVector(*M_interfaceMap.getMap(Unique)));

    for(int j=0; j<numtasks; ++j)
        epetraWorldComm->Broadcast( &numInterfaceDof[j], 1, j);

    for(int j=numtasks-1; j>0 ; --j)
    {
        numInterfaceDof[j] = numInterfaceDof[j-1];
    }
    numInterfaceDof[0]=0;
    for(int j=1; j<numtasks ; ++j)
        numInterfaceDof[j] += numInterfaceDof[j-1];

    UInt k=1;
    UInt l=0;

    M_interface = (UInt) interfaceMap.getMap(Unique)->NumGlobalElements()/nDimensions;
//UInt solidDim=M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions;
    for(l=0, ITrow=locDofMap.begin(); ITrow!=locDofMap.end() ; ++ITrow)
    {
        if(interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/)>=0)
        {
            (*M_numerationInterface)[ITrow->second /*+ dim*solidDim*/ ]=l+1+ (int)(numInterfaceDof[pid]/nDimensions)/*+ dim*localInterface*/      ;
            //                                    (*M_numerationInterfaceInt)[ITrow->second /*+ dim*solidDim*/ ]=l+1+ (int)(M_numInterfaceDof[pid]/nDimensions)/*+ dim*localInterface*/      ;
            if((int)(*M_numerationInterface)(ITrow->second )!=floor(l+1+ numInterfaceDof[pid]/nDimensions+0.2 /*+ dim*localInterface*/) )
                std::cout<<"ERROR! the numeration of the coupling map is not correct"<<std::endl;
            ++l;
        }
    }

    std::vector<int> couplingVector;
    couplingVector.reserve((int)(interfaceMap.getMap(Unique)->NumMyElements()));

    for(int dim=0; dim<nDimensions; ++dim)
    {
        for( ITrow=locDofMap.begin(); ITrow!=locDofMap.end() ; ++ITrow)
        {
            if(interfaceMap.getMap(Unique)->LID(ITrow->second)>=0)
            {
                couplingVector.push_back((*M_numerationInterface)(ITrow->second /*+ dim * solidDim*/)+ dim * M_interface );
                //couplingVector.push_back((*M_numerationInterfaceInt)[ITrow->second /*+ dim * solidDim*/]+ dim * M_interface );
            }
        }
    }// so the map for the coupling part of the matrix is just Unique

    M_interfaceMap.reset(new EpetraMap(-1, static_cast< Int> ( couplingVector.size() ), &couplingVector[0], interfaceMap.getMap(Repeated)->IndexBase()/*1*/, epetraWorldComm));
}

void BlockMatrix::applyBoundaryConditions(const Real& time)
{
    for( UInt i = 0; i < super::M_blocks.size(); ++i )
        applyBoundaryConditions( time, i );
}

void BlockMatrix::applyBoundaryConditions(const Real& time, vector_ptrtype& rhs)
{
    for( UInt i = 0; i < super::M_blocks.size(); ++i )
        applyBoundaryConditions( time, rhs, i );
}

void BlockMatrix::applyBoundaryConditions(const Real& time, vector_ptrtype& rhs, const UInt block)
{
    bcManage( *M_globalMatrix , *rhs, *super::M_FESpace[block]->mesh(), super::M_FESpace[block]->dof(), *super::M_bch[block], super::M_FESpace[block]->feBd(), 1., time);
}

void BlockMatrix::applyBoundaryConditions(const Real& time, const UInt block)
{
    bcManageMatrix( *M_globalMatrix , *super::M_FESpace[block]->mesh(), super::M_FESpace[block]->dof(), *super::M_bch[block], super::M_FESpace[block]->feBd(), 1., time);
}

void BlockMatrix::addToCoupling( const matrix_ptrtype& Mat, UInt /*position*/)
{
    if(!M_coupling->getMatrixPtr()->Filled())
        *M_coupling += *Mat;
    else
    {
        matrix_ptrtype tmp(new matrix_type(M_coupling->getMap()));
        *tmp += *M_coupling;
        *tmp += *Mat;
        M_coupling = tmp;
    }
}

} // Namespace LifeV
