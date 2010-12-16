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

#include <BlockMatrixRN.hpp>
;
namespace LifeV {


// ===================================================
//! Public Methods
// ===================================================


void BlockMatrixRN::setDataFromGetPot( const GetPot& data, const std::string& section )
{
    super::setDataFromGetPot(data, section);
    setRobinData(data, section + "/robin");
}

void BlockMatrixRN::coupler(map_shared_ptrtype& map,
                            const std::map<ID, ID>& locDofMap,
                            const vectorPtr_Type& numerationInterface,
                            const Real& timeStep)
{
    super::coupler( map,/* M_FESpace[0], M_offset[0], M_FESpace[1], M_offset[1],*/ locDofMap, numerationInterface, timeStep );
    M_robinCoupling.reset(new matrix_Type(M_coupling->getMap(), 0));
    robinCoupling( M_robinCoupling, M_alphaf, M_alphas, 7, M_FESpace[1], M_offset[1], M_FESpace[0], M_offset[0], locDofMap, numerationInterface );
    M_robinCoupling->GlobalAssemble( );
    //    M_robinCoupling->spy("RC");
}

void BlockMatrixRN::GlobalAssemble()
{
    super::GlobalAssemble();
    vector_Type rhs((*M_robinCoupling)*(*M_rhsVec));
    *M_rhsVec += rhs;
}

void BlockMatrixRN::blockAssembling()
{
    M_coupling->GlobalAssemble();
    M_globalMatrix.reset(new matrix_Type(M_coupling->getMap()));
    *M_globalMatrix += *M_coupling;
    for (UInt k=0; k<super::M_blocks.size(); ++k)
    {
        super::M_blocks[k]->GlobalAssemble();
    }

    applyRobinCoupling(M_blocks);
    M_robinPart->GlobalAssemble();

    matrixPtr_Type fluidRobinBlock(new matrix_Type(M_robinCoupling->getMap(), 1));
    *fluidRobinBlock += *M_blocks[1];
    *fluidRobinBlock += *M_robinPart;
    fluidRobinBlock->GlobalAssemble();
    M_blocks[1]->swapCrsMatrix( fluidRobinBlock->getMatrixPtr() );

    for (UInt k=0; k<M_blocks.size(); ++k)
    {
        *M_globalMatrix += *M_blocks[k];
    }
}

} // Namespace LifeV
