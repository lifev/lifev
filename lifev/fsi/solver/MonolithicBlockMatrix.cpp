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

#include <lifev/core/LifeV.hpp>

#include <lifev/fsi/solver/MonolithicBlockMatrix.hpp>

namespace LifeV
{



// ===================================================
//! Public Methods
// ===================================================


void MonolithicBlockMatrix::setDataFromGetPot (const GetPot& /*data*/, const std::string& /*section*/)
{
}

void MonolithicBlockMatrix::setupSolver (solver_Type& solver, const GetPot& data)
{
    solver.setupPreconditioner (data, "linear_system/prec"); //to avoid if we have already a prec.
}


void MonolithicBlockMatrix::GlobalAssemble()
{
    M_globalMatrix->globalAssemble();
    //    M_globalMatrix->spy("Prec");
}

void MonolithicBlockMatrix::coupler (mapPtr_Type& map,
                                     const std::map<ID, ID>& locDofMap,
                                     const vectorPtr_Type& numerationInterface,
                                     const Real& timeStep,
                                     const Real& coefficient = 1.,
                                     const Real& rescaleFactor = 1.
                                    )
{
    ASSERT (!M_coupling.get(), "coupler must not be called twice \n");
    M_coupling.reset (new matrix_Type (*map) );
    super_Type::couplingMatrix ( M_coupling,  (*M_couplingFlags) [0], super_Type::M_FESpace, super_Type::M_offset, locDofMap, numerationInterface, timeStep, 1.,  coefficient, rescaleFactor);
}


void MonolithicBlockMatrix::coupler (mapPtr_Type& /*map*/,
                                     const std::map<ID, ID>& locDofMap,
                                     const vectorPtr_Type& numerationInterface,
                                     const Real& timeStep,
                                     const Real& coefficient,
                                     const Real& rescaleFactor,
                                     UInt couplingBlock
                                    )
{
    super_Type::couplingMatrix ( M_coupling,  (*M_couplingFlags) [couplingBlock], super_Type::M_FESpace, super_Type::M_offset, locDofMap, numerationInterface, timeStep, 1., coefficient, rescaleFactor);
}


int MonolithicBlockMatrix::solveSystem ( const vector_Type& rhs, vector_Type& step, solverPtr_Type& linearSolver)
{
    return linearSolver->solveSystem (rhs, step, M_globalMatrix);
}

void MonolithicBlockMatrix::push_back_matrix ( const matrixPtr_Type& Mat, bool /*recompute*/)
{
    super_Type::M_blocks.push_back (Mat);
}


void MonolithicBlockMatrix::replace_matrix ( const matrixPtr_Type& Mat, UInt index)
{
    super_Type::M_blocks[index] = Mat;
}


void MonolithicBlockMatrix::replace_precs ( const epetraOperatorPtr_Type& /*Mat*/, UInt /*index*/)
{
    assert (false);
}


void MonolithicBlockMatrix::blockAssembling()
{
    M_coupling->globalAssemble();
    M_globalMatrix.reset (new matrix_Type (M_coupling->map() ) );
    *M_globalMatrix += *M_coupling;
    for (UInt k = 0; k < M_blocks.size(); ++k)
    {
        M_blocks[k]->globalAssemble();
        *M_globalMatrix += *M_blocks[k];
    }
}



void MonolithicBlockMatrix::applyPreconditioner (const matrixPtr_Type matrix, vectorPtr_Type& rhsFull)
{
    this->applyPreconditioner (matrix, M_globalMatrix);
    *rhsFull = (*matrix) * (*rhsFull);
}


void MonolithicBlockMatrix::applyPreconditioner ( matrixPtr_Type robinCoupling, matrixPtr_Type prec, vectorPtr_Type& rhs)
{
    applyPreconditioner (robinCoupling, prec);
    applyPreconditioner (robinCoupling, M_globalMatrix);
    (*rhs) = (*robinCoupling) * (*rhs);
}


void MonolithicBlockMatrix::applyPreconditioner ( const matrixPtr_Type prec, matrixPtr_Type& oper )
{
    matrix_Type tmpMatrix (prec->map(), 1);
    EpetraExt::MatrixMatrix::Multiply ( *prec->matrixPtr(),
                                        false,
                                        *oper->matrixPtr(),
                                        false,
                                        *tmpMatrix.matrixPtr() );
    oper->swapCrsMatrix (tmpMatrix);
}


void MonolithicBlockMatrix::createInterfaceMap ( const MapEpetra& interfaceMap , const std::map<ID, ID>& locDofMap, const UInt subdomainMaxId,  const boost::shared_ptr<Epetra_Comm> epetraWorldComm )
{
    //std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
    std::map<ID, ID>::const_iterator ITrow;

    Int numtasks = epetraWorldComm->NumProc();
    int* numInterfaceDof (new int[numtasks]);
    int pid = epetraWorldComm->MyPID();
    int numMyElements = interfaceMap.map (Unique)->NumMyElements();
    numInterfaceDof[pid] = numMyElements;
    MapEpetra subMap (*interfaceMap.map (Unique), (UInt) 0, subdomainMaxId);

    M_numerationInterface.reset (new vector_Type (subMap, Unique) );
    //should be an int vector instead of double
    //                    M_numerationInterfaceInt.reset(new Epetra_IntVector(*M_interfaceMap.getMap(Unique)));

    for (int j = 0; j < numtasks; ++j)
    {
        epetraWorldComm->Broadcast ( &numInterfaceDof[j], 1, j);
    }

    for (int j = numtasks - 1; j > 0 ; --j)
    {
        numInterfaceDof[j] = numInterfaceDof[j - 1];
    }
    numInterfaceDof[0] = 0;
    for (int j = 1; j < numtasks ; ++j)
    {
        numInterfaceDof[j] += numInterfaceDof[j - 1];
    }

    UInt l = 0;

    M_interface = (UInt) interfaceMap.map (Unique)->NumGlobalElements() / nDimensions;
    //UInt solidDim=M_dFESpace->map().map(Unique)->NumGlobalElements()/nDimensions;
    for (l = 0, ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
    {
        if (interfaceMap.map (Unique)->LID (ITrow->second /*+ dim*solidDim*/) >= 0)
        {
            (*M_numerationInterface) [ITrow->second /*+ dim*solidDim*/ ] = l + (int) (numInterfaceDof[pid] / nDimensions) /*+ dim*localInterface*/      ;
            //                                    (*M_numerationInterfaceInt)[ITrow->second /*+ dim*solidDim*/ ]=l+1+ (int)(M_numInterfaceDof[pid]/nDimensions)/*+ dim*localInterface*/      ;
            if ( (int) (*M_numerationInterface) (ITrow->second ) != floor (l + numInterfaceDof[pid] / nDimensions + 0.2 /*+ dim*localInterface*/) )
            {
                std::cout << "ERROR! the numeration of the coupling map is not correct" << std::endl;
            }
            ++l;
        }
    }

    std::vector<int> couplingVector;
    couplingVector.reserve ( (int) (interfaceMap.map (Unique)->NumMyElements() ) );

    for (UInt dim = 0; dim < nDimensions; ++dim)
    {
        for ( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
        {
            if (interfaceMap.map (Unique)->LID (ITrow->second) >= 0)
            {
                couplingVector.push_back ( (*M_numerationInterface) (ITrow->second /*+ dim * solidDim*/) + dim * M_interface );
                //couplingVector.push_back((*M_numerationInterfaceInt)[ITrow->second /*+ dim * solidDim*/]+ dim * M_interface );
            }
        }
    }// so the map for the coupling part of the matrix is just Unique

    M_interfaceMap.reset (new MapEpetra (-1, static_cast< Int> ( couplingVector.size() ), &couplingVector[0], epetraWorldComm) );
}

void MonolithicBlockMatrix::applyBoundaryConditions (const Real& time)
{
    for ( UInt i = 0; i < super_Type::M_blocks.size(); ++i )
    {
        applyBoundaryConditions ( time, i );
    }
}

void MonolithicBlockMatrix::applyBoundaryConditions (const Real& time, vectorPtr_Type& rhs)
{
    for ( UInt i = 0; i < super_Type::M_blocks.size(); ++i )
    {
        applyBoundaryConditions ( time, rhs, i );
    }
}

void MonolithicBlockMatrix::applyBoundaryConditions (const Real& time, vectorPtr_Type& rhs, const UInt block)
{

    std::cout << "Block: " <<  block << std::endl;
    ( *super_Type::M_bch[block] ).showMe();

    for ( ID i = 0; i < ( *super_Type::M_bch[block] ).size(); ++i )
    {
        if ( ( ( *super_Type::M_bch[block] )[i]).isDataAVector() )
        {
            std::cout << "Size:" <<  ( ( *super_Type::M_bch[block] )[i]).pointerToBCVector()->rhsVector().map().mapSize() << std::endl;
        }
    }

    bcManage ( *M_globalMatrix , *rhs, *super_Type::M_FESpace[block]->mesh(), super_Type::M_FESpace[block]->dof(), *super_Type::M_bch[block], super_Type::M_FESpace[block]->feBd(), 1., time);
}

void MonolithicBlockMatrix::applyBoundaryConditions (const Real& time, const UInt block)
{

    bcManageMatrix ( *M_globalMatrix , *super_Type::M_FESpace[block]->mesh(), super_Type::M_FESpace[block]->dof(), *super_Type::M_bch[block], super_Type::M_FESpace[block]->feBd(), 1., time);
}

void MonolithicBlockMatrix::addToCoupling ( const matrixPtr_Type& Mat, UInt /*position*/)
{
    if (!M_coupling->matrixPtr()->Filled() )
    {
        *M_coupling += *Mat;
    }
    else
    {
        matrixPtr_Type tmp (new matrix_Type (M_coupling->map() ) );
        *tmp += *M_coupling;
        *tmp += *Mat;
        M_coupling = tmp;
    }
}


void MonolithicBlockMatrix::addToCoupling ( const Real& entry , UInt row, UInt col, UInt /*position*/ )
{
    if (!M_coupling->matrixPtr()->Filled() )
    {
        M_coupling->setCoefficient (row, col, entry);
    }
    else
    {
        matrixPtr_Type tmp (new matrix_Type (M_coupling->map() ) );
        *tmp += *M_coupling;
        tmp->setCoefficient (row, col, entry);
        M_coupling = tmp;
    }
}


} // Namespace LifeV
