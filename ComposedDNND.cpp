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

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 02 Jul 2010
 */

#include <ComposedDNND.hpp>

namespace LifeV {

void ComposedDNND::coupler(map_shared_ptrtype& map,
                     const std::map<ID, ID>& locDofMap,
                           const vector_ptrtype& numerationInterface,
                           const Real& timeStep)
{
    UInt totalDofs=map->getMap(Unique)->NumGlobalElements()+1;
    UInt fluidSolid=M_offset[0]+1+M_FESpace[0]->map().getMap(Unique)->NumGlobalElements();

    for(ID k=0; k<2; ++k)
    {
        M_blocks[k]->GlobalAssemble();
        matrix_ptrtype block(new matrix_type(*M_blocks[k]));
        M_blocks.push_back(block);
        M_bch.push_back(M_bch[k]);
        M_FESpace.push_back(M_FESpace[k]);
        M_offset.push_back(M_offset[k]);
        M_recompute[2+k]=(M_recompute[k]);
        M_blockReordering->push_back((*M_blockReordering)[k]);
    }



    matrix_ptrtype coupling(new matrix_type(*map));
    UInt one(1.);

    coupling.reset(new matrix_type(*map, 0));
    coupling->insertValueDiagonal(one, M_offset[fluid]+1, M_offset[solid]+1 );
    coupling->insertValueDiagonal(one, fluidSolid, totalDofs);
    couplingMatrix(coupling, (*M_couplingFlags)[0]/*8*/, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2.);
    M_coupling.push_back(coupling);

    coupling.reset(new matrix_type(*map, 0));
    coupling->insertValueDiagonal( one, M_FESpace[solid]->map() , M_offset[solid] );
    coupling->insertValueDiagonal(one, fluidSolid, totalDofs);
    couplingMatrix(coupling, (*M_couplingFlags)[1]/*4*/, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 4./**/);
    couplingMatrix(coupling, (*M_couplingFlags)[2]/*2*/, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2.);
    M_coupling.push_back(coupling);

    coupling.reset( new matrix_type( *map, 0 ) );
    coupling->insertValueDiagonal( one,  M_offset[fluid]+1, M_offset[solid]+1 );
    coupling->insertValueDiagonal( -1, fluidSolid, totalDofs );
    couplingMatrix( coupling, (*M_couplingFlags)[3]/*8*/,  M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2. );
    couplingMatrix( coupling, (*M_couplingFlags)[4]/*1*/,  M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 4./**/ );
    M_coupling.push_back( coupling );

    coupling.reset(new matrix_type( *map, 0 ));
    coupling->insertValueDiagonal( one, M_FESpace[solid]->map(), M_offset[solid] );
    coupling->insertValueDiagonal( one, fluidSolid, totalDofs );
    couplingMatrix( coupling, (*M_couplingFlags)[5]/*2*/, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2. );
    M_coupling.push_back( coupling );

    M_prec.resize(M_blocks.size());
    //M_blockPrecs.reset(new ComposedPreconditioner<ComposedPreconditioner<Ifpack_Preconditioner> >(&M_blocks[0]->getMatrixPtr()->Comm()));
}

} // Namespace LifeV
