/* -*- mode: c++ -*- */
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
    @brief A short description of the file content

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 08 Jun 2010

    A more detailed description of the file (if necessary)
 */

#include <life/lifecore/life.hpp>

#include <BlockMatrix.hpp>

namespace LifeV
{



// ===================================================
//! Public Methods
// ===================================================


void BlockMatrix::setDataFromGetPot(const GetPot& data, const std::string& section)
{
}


void BlockMatrix::GlobalAssemble()
{
    M_globalMatrix->GlobalAssemble();
    //    M_globalMatrix->spy("Prec");
}

void BlockMatrix::coupler(mapPtr_Type& map,
                          const std::map<ID, ID>& locDofMap,
                          const vectorPtr_Type& numerationInterface,
                          const Real& timeStep
                         )
{
    ASSERT(!M_coupling.get(), "coupler must not be called twice \n");
    M_coupling.reset(new matrix_Type(*map));
    super_Type::couplingMatrix( M_coupling,  M_couplingFlag, super_Type::M_FESpace, super_Type::M_offset, locDofMap, numerationInterface, timeStep);
}

void BlockMatrix::coupler(mapPtr_Type& map,
                          const std::map<ID, ID>& locDofMap,
                          const vectorPtr_Type& numerationInterface,
                          const Real& timeStep,
                          UInt /*flag1*/
                         )
{
    super_Type::couplingMatrix( M_coupling,  M_couplingFlag, super_Type::M_FESpace, super_Type::M_offset, locDofMap, numerationInterface, timeStep);
}


int BlockMatrix::solveSystem( const vector_Type& rhs, vector_Type& step, solverPtr_Type& linearSolver)
{
    return linearSolver->solveSystem(rhs, step, M_globalMatrix);
}

void BlockMatrix::push_back_matrix( const matrixPtr_Type& Mat, bool /*recompute*/)
{
    super_Type::M_blocks.push_back(Mat);
}


void BlockMatrix::replace_matrix( const matrixPtr_Type& Mat, UInt index)
{
    super_Type::M_blocks[index] = Mat;
}


void BlockMatrix::replace_precs( const epetraOperatorPtr_Type& Mat, UInt index)
{
    assert(false);
}


void BlockMatrix::blockAssembling()
{
    M_coupling->GlobalAssemble();
    M_globalMatrix.reset(new matrix_Type(M_coupling->getMap()));
    *M_globalMatrix += *M_coupling;
    for (UInt k=0; k<M_blocks.size(); ++k)
    {
        M_blocks[k]->GlobalAssemble();
        *M_globalMatrix += *M_blocks[k];
    }
}



void BlockMatrix::applyPreconditioner(const matrixPtr_Type matrix, vectorPtr_Type& rhsFull)
{
    this->applyPreconditioner(matrix, M_globalMatrix);
    *rhsFull = (*matrix)*(*rhsFull);
}


void BlockMatrix::applyPreconditioner( matrixPtr_Type robinCoupling, matrixPtr_Type prec, vectorPtr_Type& rhs)
{
    applyPreconditioner(robinCoupling, prec);
    applyPreconditioner(robinCoupling, M_globalMatrix);
    (*rhs) = (*robinCoupling)*(*rhs);
}


void BlockMatrix::applyPreconditioner( const matrixPtr_Type prec, matrixPtr_Type& oper )
{
    matrix_Type tmpMatrix(prec->getMap(), 1);
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

    M_numerationInterface.reset(new vector_Type(subMap,Unique));
    //should be an int vector instead of double
    //                    M_numerationInterfaceInt.reset(new Epetra_IntVector(*M_interfaceMap.getMap(Unique)));

    for (int j=0; j<numtasks; ++j)
        epetraWorldComm->Broadcast( &numInterfaceDof[j], 1, j);

    for (int j=numtasks-1; j>0 ; --j)
    {
        numInterfaceDof[j] = numInterfaceDof[j-1];
    }
    numInterfaceDof[0]=0;
    for (int j=1; j<numtasks ; ++j)
        numInterfaceDof[j] += numInterfaceDof[j-1];

    UInt k=1;
    UInt l=0;

    M_interface = (UInt) interfaceMap.getMap(Unique)->NumGlobalElements()/nDimensions;
//UInt solidDim=M_dFESpace->map().getMap(Unique)->NumGlobalElements()/nDimensions;
    for (l=0, ITrow=locDofMap.begin(); ITrow!=locDofMap.end() ; ++ITrow)
    {
        if (interfaceMap.getMap(Unique)->LID(ITrow->second /*+ dim*solidDim*/)>=0)
        {
            (*M_numerationInterface)[ITrow->second /*+ dim*solidDim*/ ]=l+1+ (int)(numInterfaceDof[pid]/nDimensions)/*+ dim*localInterface*/      ;
            //                                    (*M_numerationInterfaceInt)[ITrow->second /*+ dim*solidDim*/ ]=l+1+ (int)(M_numInterfaceDof[pid]/nDimensions)/*+ dim*localInterface*/      ;
            if ((int)(*M_numerationInterface)(ITrow->second )!=floor(l+1+ numInterfaceDof[pid]/nDimensions+0.2 /*+ dim*localInterface*/) )
                std::cout<<"ERROR! the numeration of the coupling map is not correct"<<std::endl;
            ++l;
        }
    }

    std::vector<int> couplingVector;
    couplingVector.reserve((int)(interfaceMap.getMap(Unique)->NumMyElements()));

    for (int dim=0; dim<nDimensions; ++dim)
    {
        for ( ITrow=locDofMap.begin(); ITrow!=locDofMap.end() ; ++ITrow)
        {
            if (interfaceMap.getMap(Unique)->LID(ITrow->second)>=0)
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
    for ( UInt i = 0; i < super_Type::M_blocks.size(); ++i )
        applyBoundaryConditions( time, i );
}

void BlockMatrix::applyBoundaryConditions(const Real& time, vectorPtr_Type& rhs)
{
    for ( UInt i = 0; i < super_Type::M_blocks.size(); ++i )
        applyBoundaryConditions( time, rhs, i );
}

void BlockMatrix::applyBoundaryConditions(const Real& time, vectorPtr_Type& rhs, const UInt block)
{
    bcManage( *M_globalMatrix , *rhs, *super_Type::M_FESpace[block]->mesh(), super_Type::M_FESpace[block]->dof(), *super_Type::M_bch[block], super_Type::M_FESpace[block]->feBd(), 1., time);
}

void BlockMatrix::applyBoundaryConditions(const Real& time, const UInt block)
{
    bcManageMatrix( *M_globalMatrix , *super_Type::M_FESpace[block]->mesh(), super_Type::M_FESpace[block]->dof(), *super_Type::M_bch[block], super_Type::M_FESpace[block]->feBd(), 1., time);
}

void BlockMatrix::addToCoupling( const matrixPtr_Type& Mat, UInt /*position*/)
{
    if (!M_coupling->getMatrixPtr()->Filled())
        *M_coupling += *Mat;
    else
    {
        matrixPtr_Type tmp(new matrix_Type(M_coupling->getMap()));
        *tmp += *M_coupling;
        *tmp += *Mat;
        M_coupling = tmp;
    }
}

} // Namespace LifeV
