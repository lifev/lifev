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

#include <ComposedDN.hpp>

namespace LifeV
{


// ===================================================
//! Public Methods
// ===================================================

void ComposedDN::setDataFromGetPot( const GetPot& dataFile,
                                      const std::string& section )
{
    M_blockPrecs->setDataFromGetPot( dataFile, section );
}


int ComposedDN::solveSystem( const vector_Type& rhs, vector_Type& step, solver_ptrtype& linearSolver )
{
    assert(M_blockPrecs.get());

    if (!set())
    {
        for (UInt k=0; k < M_blocks.size(); ++k)
            push_back_precs(M_blocks[(*M_blockReordering)[k]]);
    }
    else
    {
        for (UInt k=0; k < M_blocks.size(); ++k)
        {
            if (M_recompute[(*M_blockReordering)[k]])
            {
                linearSolver->displayer()->leaderPrint("  M-  Computing preconditioner factor:           ", k, "\n");
                replace_precs(M_blocks[(*M_blockReordering)[k]], k);
            }
            else
                linearSolver->displayer()->leaderPrint("  M-  Reusing preconditioner factor:           ", k, "\n");
        }
    }
    return linearSolver->solveSystem(rhs, step, boost::static_pointer_cast<EpetraPreconditioner>(M_blockPrecs));
}


void ComposedDN::coupler(map_shared_ptrtype& map,
                         const std::map<ID, ID>& locDofMap,
                         const vectorPtr_Type& numerationInterface,
                         const Real& timeStep)
{
    UInt totalDofs( map->getMap(Unique)->NumGlobalElements() );
    UInt solidAndFluid(M_offset[solid]+1+M_FESpace[solid]->map().getMap(Unique)->NumGlobalElements());

    matrixPtr_Type coupling(new matrix_Type(*map));
    couplingMatrix( coupling,  (*M_couplingFlags)[solid], M_FESpace, M_offset, locDofMap, numerationInterface, timeStep);
    coupling->insertValueDiagonal( 1., M_offset[fluid]+1, M_offset[solid]+1 );
    coupling->insertValueDiagonal( 1., solidAndFluid, totalDofs+1 );
    M_coupling.push_back(coupling);

    coupling.reset(new matrix_Type(*map));
    couplingMatrix( coupling,  (*M_couplingFlags)[fluid], M_FESpace, M_offset, locDofMap, numerationInterface, timeStep);
    coupling->insertValueDiagonal( 1. , M_offset[solid], solidAndFluid );
    coupling->insertValueDiagonal( 1. , solidAndFluid + nDimensions*numerationInterface->getMap().getMap(Unique)->NumGlobalElements(), totalDofs +1 );
    M_coupling.push_back(coupling);

}

void ComposedDN::push_back_precs( matrixPtr_Type& Mat )
{
    M_blockPrecs->push_back(Mat);
}


// ===================================================
//! Protected Methods
// ===================================================


void ComposedDN::replace_precs(  matrixPtr_Type& Mat, UInt position )
{
    M_blockPrecs->replace(Mat, position);
}


} // Namespace LifeV
