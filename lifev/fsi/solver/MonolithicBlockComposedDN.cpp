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

#include <lifev/fsi/solver/MonolithicBlockComposedDN.hpp>

namespace LifeV
{


// ===================================================
//! Public Methods
// ===================================================

void MonolithicBlockComposedDN::setDataFromGetPot( const GetPot& dataFile,
                                      const std::string& section )
{
    M_blockPrecs->setDataFromGetPot( dataFile, section );
}


int MonolithicBlockComposedDN::solveSystem( const vector_Type& rhs, vector_Type& step, solverPtr_Type& linearSolver )
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
                linearSolver->displayer()->leaderPrint("  M-  Computing preconditioner factor:         ", k, "\n");
                replace_precs(M_blocks[(*M_blockReordering)[k]], k);
            }
            else
                linearSolver->displayer()->leaderPrint("  M-  Reusing preconditioner factor:           ", k, "\n");
        }
    }
    return linearSolver->solveSystem(rhs, step, boost::static_pointer_cast<Preconditioner>(M_blockPrecs));
}


void MonolithicBlockComposedDN::coupler(mapPtr_Type& map,
                                        const std::map<ID, ID>& locDofMap,
                                        const vectorPtr_Type& numerationInterface,
                                        const Real& timeStep,
                                        const Real& coefficient,
                                        const Real& rescaleFactor)
{
    UInt totalDofs( map->map(Unique)->NumGlobalElements() );
    UInt solidAndFluid(M_offset[solid]+M_FESpace[solid]->map().map(Unique)->NumGlobalElements());

    matrixPtr_Type coupling(new matrix_Type(*map));
    couplingMatrix( coupling,  (*M_couplingFlags)[solid], M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 1., coefficient, rescaleFactor);
    coupling->insertValueDiagonal( 1., M_offset[fluid], M_offset[solid] );
    coupling->insertValueDiagonal( 1., solidAndFluid, totalDofs );
    M_coupling.push_back(coupling);

    coupling.reset(new matrix_Type(*map));
    couplingMatrix( coupling,  (*M_couplingFlags)[fluid], M_FESpace, M_offset, locDofMap, numerationInterface, timeStep, 1., coefficient, rescaleFactor);
    coupling->insertValueDiagonal( 1. , M_offset[solid], solidAndFluid );
    coupling->insertValueDiagonal( 1. , solidAndFluid + nDimensions*numerationInterface->map().map(Unique)->NumGlobalElements(), totalDofs );
    M_coupling.push_back(coupling);

}

void MonolithicBlockComposedDN::push_back_precs( matrixPtr_Type& Mat )
{
    M_blockPrecs->push_back(Mat);
}


// ===================================================
//! Protected Methods
// ===================================================


void MonolithicBlockComposedDN::replace_precs(  matrixPtr_Type& Mat, UInt position )
{
    M_blockPrecs->replace(Mat, position);
}


} // Namespace LifeV
