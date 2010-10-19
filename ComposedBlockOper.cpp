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
    @file ComposedBlockOper.cpp

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 08 Jun 2010

 */

#include <ComposedBlockOper.hpp>

namespace LifeV {

void ComposedBlockOper::blockAssembling()
{
    for(UInt k=0; k<M_blocks.size(); ++k)
    {
        blockAssembling(k);
    }
//     M_coupling[0]->spy("C1");
//     M_coupling[1]->spy("C2");
}

void ComposedBlockOper::blockAssembling(const UInt k)
{
    if(M_blocks[k]->getMatrixPtr()->Filled())
    {
        matrix_ptrtype tmp(new matrix_type(M_blocks[0]->getMap(), 1));
        *tmp += *M_blocks[k];
        M_blocks[k]=tmp;
    }
    M_coupling[k]->GlobalAssemble();
    *M_blocks[k] += *M_coupling[k];
}

void ComposedBlockOper::GlobalAssemble()
{
    for(UInt k=0; k<M_blocks.size(); ++k)
    {
        M_blocks[k]->GlobalAssemble();
    }
//             M_blocks[0]->spy("first");
//             M_blocks[1]->spy("second");
//             M_blocks[2]->spy("third");
//        M_blocks[3]->spy("fourth");
}

void ComposedBlockOper::push_back_matrix(const matrix_ptrtype& Mat, const  bool recompute)
{
    M_blocks.push_back(Mat);
    M_recompute[M_blocks.size()-1] = recompute;
}

void ComposedBlockOper::replace_matrix( const matrix_ptrtype& Mat, UInt position )
{
    M_blocks[position]=Mat;
}

void ComposedBlockOper::replace_coupling( const matrix_ptrtype& Mat, UInt position )
{
    M_coupling[position]=Mat;
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

void ComposedBlockOper::addToCoupling( const matrix_ptrtype& Mat, UInt position)
{
     Mat->GlobalAssemble();
    *M_coupling[position] += *Mat;
}

void ComposedBlockOper::push_back_oper( ComposedBlockOper& Oper)
{
    super::push_back_oper( Oper );
    M_coupling.insert(M_coupling.end(), Oper.getCouplingVector().begin(), Oper.getCouplingVector().end());
}

void ComposedBlockOper::coupler( map_shared_ptrtype& map,
                                 const std::map<ID, ID>& locDofMap,
                                 const vector_ptrtype& numerationInterface,
                                 const Real& timeStep,
                                 UInt couplingBlock )
{
    matrix_ptrtype coupling(new matrix_type(*map));
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

void
ComposedBlockOper::push_back_coupling( matrix_ptrtype& coupling )
{
    M_coupling.push_back(coupling);
}

} // Namespace LifeV
