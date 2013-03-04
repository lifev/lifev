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

#include <lifev/fsi/solver/MonolithicBlockComposed.hpp>

namespace LifeV
{


// ===================================================
//! Public Methods
// ===================================================

void MonolithicBlockComposed::GlobalAssemble()
{
    for (UInt k = 0; k < M_blocks.size(); ++k)
    {
        M_blocks[k]->globalAssemble();
    }
    //     M_blocks[0]->spy("first");
    //     M_blocks[1]->spy("second");
    //     M_blocks[2]->spy("third");
    //        M_blocks[3]->spy("fourth");
}

void MonolithicBlockComposed::blockAssembling()
{
    for (UInt k = 0; k < M_blocks.size(); ++k)
    {
        blockAssembling (k);
    }
    //     M_coupling[0]->spy("C1");
    //     M_coupling[1]->spy("C2");
    //     M_coupling[2]->spy("C1");
    //     M_coupling[3]->spy("C2");
}



void MonolithicBlockComposed::coupler ( mapPtr_Type& map,
                                        const std::map<ID, ID>& locDofMap,
                                        const vectorPtr_Type& numerationInterface,
                                        const Real& timeStep,
                                        const Real& coefficient,
                                        const Real& rescaleFactor,
                                        UInt couplingBlock)
{
    matrixPtr_Type coupling (new matrix_Type (*map) );
    couplingMatrix ( coupling,  (*M_couplingFlags) [couplingBlock], M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 1., coefficient, rescaleFactor);
    UInt totalDofs ( map->map (Unique)->NumGlobalElements() );

    coupling->insertValueDiagonal ( 1., 0 , M_offset[couplingBlock] );
    coupling->insertValueDiagonal ( 1., M_offset[couplingBlock] + M_FESpace[couplingBlock]->map().map (Unique)->NumGlobalElements(), totalDofs );

    if (couplingBlock != M_coupling.size() + 1)
    {
        M_coupling.insert (M_coupling.begin() + couplingBlock, coupling);
    }
    else
    {
        M_coupling.push_back (coupling);
    }
}

void MonolithicBlockComposed::push_back_matrix (const matrixPtr_Type& Mat, const  bool recompute)
{
    M_blocks.push_back (Mat);
    M_recompute[M_blocks.size() - 1] = recompute;
}

void MonolithicBlockComposed::push_back_oper ( MonolithicBlockComposed& Oper)
{
    super_Type::push_back_oper ( Oper );
    M_coupling.insert (M_coupling.end(), Oper.couplingVector().begin(), Oper.couplingVector().end() );
}

void
MonolithicBlockComposed::push_back_coupling ( matrixPtr_Type& coupling )
{
    M_coupling.push_back (coupling);
}

void MonolithicBlockComposed::replace_matrix ( const matrixPtr_Type& Mat, UInt position )
{
    M_blocks[position] = Mat;
}

void MonolithicBlockComposed::replace_coupling ( const matrixPtr_Type& Mat, UInt position )
{
    M_coupling[position] = Mat;
}

void MonolithicBlockComposed::addToCoupling ( const matrixPtr_Type& Mat, UInt position)
{
    Mat->globalAssemble();
    *M_coupling[position] += *Mat;
}

void MonolithicBlockComposed::addToCoupling ( const Real& entry , UInt row, UInt col, UInt position )
{
    if (!M_coupling[position]->matrixPtr()->Filled() )
    {
        M_coupling[position]->setCoefficient (row, col, entry);
    }
    else
    {
        matrixPtr_Type tmp (new matrix_Type (M_coupling[position]->map() ) );
        *tmp += *M_coupling[position];
        tmp->setCoefficient (row, col, entry);
        M_coupling[position] = tmp;
    }
}

// ===================================================
//! Protected Methods
// ===================================================


void MonolithicBlockComposed::blockAssembling (const UInt k)
{
    if (M_blocks[k]->matrixPtr()->Filled() )
    {
        matrixPtr_Type tmp (new matrix_Type (M_blocks[0]->map(), 1) );
        *tmp += *M_blocks[k];
        M_blocks[k] = tmp;
    }
    M_coupling[k]->globalAssemble();
    *M_blocks[k] += *M_coupling[k];
}

void MonolithicBlockComposed::swap (const UInt i, const UInt j)
{
    super_Type::swap (M_blocks[i], M_blocks[j]);
    super_Type::swap (M_bch[i], M_bch[j]);
    super_Type::swap (M_FESpace[i], M_FESpace[j]);
    super_Type::swap (M_coupling[i], M_coupling[j]);

    bool tmpRecompute = this->M_recompute[i];
    this->M_recompute[i] = this->M_recompute[j];
    this->M_recompute[j] = tmpRecompute;

    UInt tmpOffset = this->M_offset[i];
    this->M_offset[i] = this->M_offset[j];
    this->M_offset[j] = tmpOffset;
}

const UInt MonolithicBlockComposed::whereIsBlock ( UInt position ) const
{
    for (UInt i = 0; i < M_blockReordering->size(); i++)
    {
        if ( (*M_blockReordering) [i] == position)
        {
            return i;
        }
    }
    ERROR_MSG ("requested a block that does not exist in MonolithicBlockComposed.cpp");
}

} // Namespace LifeV
