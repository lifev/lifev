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

#include <EpetraExt_MatrixMatrix.h>

#include <MonolithicBlockMatrixRN.hpp>
;
namespace LifeV {


// ===================================================
//! Public Methods
// ===================================================


void MonolithicBlockMatrixRN::setDataFromGetPot( const GetPot& data, const std::string& section )
{
    super_Type::setDataFromGetPot(data, section);
    setRobinData(data, section + "/robin");
}

void MonolithicBlockMatrixRN::coupler(mapPtr_Type& map,
                            const std::map<ID, ID>& locDofMap,
                            const vectorPtr_Type& numerationInterface,
                            const Real& timeStep)
{
    super_Type::coupler( map,/* M_FESpace[0], M_offset[0], M_FESpace[1], M_offset[1],*/ locDofMap, numerationInterface, timeStep );
    M_robinCoupling.reset(new matrix_Type(M_coupling->map(), 0));
    robinCoupling( M_robinCoupling, M_alphaf, M_alphas, 7, M_FESpace[1], M_offset[1], M_FESpace[0], M_offset[0], locDofMap, numerationInterface );
    M_robinCoupling->globalAssemble( );
    //    M_robinCoupling->spy("RC");
}

void MonolithicBlockMatrixRN::GlobalAssemble()
{
    super_Type::GlobalAssemble();
    vector_Type rhs((*M_robinCoupling)*(*M_rhsVec));
    *M_rhsVec += rhs;
}

void MonolithicBlockMatrixRN::blockAssembling()
{
    M_coupling->globalAssemble();
    M_globalMatrix.reset(new matrix_Type(M_coupling->map()));
    *M_globalMatrix += *M_coupling;
    for (UInt k=0; k<super_Type::M_blocks.size(); ++k)
    {
        super_Type::M_blocks[k]->globalAssemble();
    }

    applyRobinCoupling(M_blocks);
    M_robinPart->globalAssemble();

    matrixPtr_Type fluidRobinBlock(new matrix_Type(M_robinCoupling->map(), 1));
    *fluidRobinBlock += *M_blocks[1];
    *fluidRobinBlock += *M_robinPart;
    fluidRobinBlock->globalAssemble();
    M_blocks[1]->swapCrsMatrix( fluidRobinBlock->matrixPtr() );

    for (UInt k=0; k<M_blocks.size(); ++k)
    {
        *M_globalMatrix += *M_blocks[k];
    }
}

} // Namespace LifeV
