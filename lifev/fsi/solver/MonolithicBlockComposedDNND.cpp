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

#include <lifev/fsi/solver/MonolithicBlockComposedDNND.hpp>

namespace LifeV
{


// ===================================================
//! Public Methods
// ===================================================

void MonolithicBlockComposedDNND::coupler(mapPtr_Type& map,
                                          const std::map<ID, ID>& locDofMap,
                                          const vectorPtr_Type& numerationInterface,
                                          const Real& timeStep,
                                          const Real& coefficient,
                                          const Real& rescaleFactor)
{
    UInt totalDofs=map->map(Unique)->NumGlobalElements();
    UInt fluidSolid=M_offset[0]+M_FESpace[0]->map().map(Unique)->NumGlobalElements();

    for (ID k=0; k<2; ++k)
    {
        M_blocks[k]->globalAssemble();
        matrixPtr_Type block(new matrix_Type(*M_blocks[k]));
        M_blocks.push_back(block);
        M_bch.push_back(M_bch[k]);
        M_FESpace.push_back(M_FESpace[k]);
        M_offset.push_back(M_offset[k]);
        M_recompute[2+k]=(M_recompute[k]);
        M_blockReordering->push_back((*M_blockReordering)[k]);
    }



    matrixPtr_Type coupling(new matrix_Type(*map));
    UInt one(1.);

    coupling.reset(new matrix_Type(*map, 0));
    coupling->insertValueDiagonal(one, M_offset[fluid], M_offset[solid] );
    coupling->insertValueDiagonal(one, fluidSolid, totalDofs);
    couplingMatrix(coupling, (*M_couplingFlags)[0]/*8*/, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2., coefficient, rescaleFactor);
    M_coupling.push_back(coupling);

    coupling.reset(new matrix_Type(*map, 0));
    coupling->insertValueDiagonal( one, M_FESpace[solid]->map() , M_offset[solid] );
    coupling->insertValueDiagonal(one, fluidSolid, totalDofs);
    couplingMatrix(coupling, (*M_couplingFlags)[1]/*4*/, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 4., coefficient, rescaleFactor/**/);
    couplingMatrix(coupling, (*M_couplingFlags)[2]/*2*/, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2., coefficient, rescaleFactor);
    M_coupling.push_back(coupling);

    coupling.reset( new matrix_Type( *map, 0 ) );
    coupling->insertValueDiagonal( one,  M_offset[fluid], M_offset[solid] );
    coupling->insertValueDiagonal( -1, fluidSolid, totalDofs );
    couplingMatrix( coupling, (*M_couplingFlags)[3]/*8*/,  M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2., coefficient, rescaleFactor );
    couplingMatrix( coupling, (*M_couplingFlags)[4]/*1*/,  M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 4., coefficient, rescaleFactor/**/ );
    M_coupling.push_back( coupling );

    coupling.reset(new matrix_Type( *map, 0 ));
    coupling->insertValueDiagonal( one, M_FESpace[solid]->map(), M_offset[solid] );
    coupling->insertValueDiagonal( one, fluidSolid, totalDofs );
    couplingMatrix( coupling, (*M_couplingFlags)[5]/*2*/, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 2., coefficient, rescaleFactor );
    M_coupling.push_back( coupling );

    M_prec.resize(M_blocks.size());
    //M_blockPrecs.reset(new ComposedOperator<ComposedOperator<Ifpack_Preconditioner> >(&M_blocks[0]->getMatrixPtr()->Comm()));
}

} // Namespace LifeV
