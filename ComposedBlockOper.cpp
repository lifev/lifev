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

#include <life/lifecore/life.hpp>

#include <ComposedBlockOper.hpp>

namespace LifeV
{


// ===================================================
//! Public Methods
// ===================================================

void ComposedBlockOper::GlobalAssemble()
{
    for (UInt k=0; k<M_blocks.size(); ++k)
    {
        M_blocks[k]->GlobalAssemble();
    }
//             M_blocks[0]->spy("first");
//             M_blocks[1]->spy("second");
//             M_blocks[2]->spy("third");
//        M_blocks[3]->spy("fourth");
}

void ComposedBlockOper::blockAssembling()
{
    for(UInt k=0; k<M_blocks.size(); ++k)
    {
        blockAssembling(k);
    }
//     M_coupling[0]->spy("C1");
//     M_coupling[1]->spy("C2");
}


void ComposedBlockOper::coupler( map_shared_ptrtype& map,
                                 const std::map<ID, ID>& locDofMap,
                                 const vector_ptrtype& numerationInterface,
                                 const Real& timeStep,
                                 UInt couplingBlock )
{
    matrixPtr_Type coupling(new matrix_Type(*map));
    couplingMatrix( coupling,  (*M_couplingFlags)[couplingBlock], M_FESpace, M_offset, locDofMap, numerationInterface, timeStep);
    UInt totalDofs( map->getMap(Unique)->NumGlobalElements() );

    coupling->insertValueDiagonal( 1., 1 , M_offset[couplingBlock]+1 );
    coupling->insertValueDiagonal( 1., M_offset[couplingBlock]+M_FESpace[couplingBlock]->map().getMap(Unique)->NumGlobalElements()+1, totalDofs+1 );

    if(couplingBlock != M_coupling.size()+1)
    {
        M_coupling.insert(M_coupling.begin()+couplingBlock, coupling);
    }
    else
        M_coupling.push_back(coupling);
}

void ComposedBlockOper::push_back_matrix(const matrixPtr_Type& Mat, const  bool recompute)
{
    M_blocks.push_back(Mat);
    M_recompute[M_blocks.size()-1] = recompute;
}

void ComposedBlockOper::push_back_oper( ComposedBlockOper& Oper)
{
    super::push_back_oper( Oper );
    M_coupling.insert(M_coupling.end(), Oper.getCouplingVector().begin(), Oper.getCouplingVector().end());
}

void
ComposedBlockOper::push_back_coupling( matrixPtr_Type& coupling )
{
    M_coupling.push_back(coupling);
}

void ComposedBlockOper::replace_matrix( const matrixPtr_Type& Mat, UInt position )
{
    M_blocks[position]=Mat;
}

void ComposedBlockOper::replace_coupling( const matrixPtr_Type& Mat, UInt position )
{
    M_coupling[position]=Mat;
}

void ComposedBlockOper::addToCoupling( const matrixPtr_Type& Mat, UInt position)
{
     Mat->GlobalAssemble();
    *M_coupling[position] += *Mat;
}



// ===================================================
//! Protected Methods
// ===================================================


void ComposedBlockOper::blockAssembling(const UInt k)
{
    if(M_blocks[k]->getMatrixPtr()->Filled())
    {
        matrixPtr_Type tmp(new matrix_Type(M_blocks[0]->getMap(), 1));
        *tmp += *M_blocks[k];
        M_blocks[k]=tmp;
    }
    M_coupling[k]->GlobalAssemble();
    *M_blocks[k] += *M_coupling[k];
}

void ComposedBlockOper::swap(const UInt i, const UInt j)
{
    super::swap(M_blocks[i], M_blocks[j]);
    super::swap(M_bch[i], M_bch[j]);
    super::swap(M_FESpace[i], M_FESpace[j]);
    super::swap(M_coupling[i], M_coupling[j]);

    bool tmpRecompute = this->M_recompute[i];
    this->M_recompute[i] = this->M_recompute[j];
    this->M_recompute[j] = tmpRecompute;

    UInt tmpOffset = this->M_offset[i];
    this->M_offset[i] = this->M_offset[j];
    this->M_offset[j] = tmpOffset;
}


} // Namespace LifeV
