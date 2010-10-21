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
  @file ComposedDN.cpp

  @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
  @date 08 Jun 2010
*/

#include <ComposedDN.hpp>

namespace LifeV {


void ComposedDN::setDataFromGetPot( const GetPot& dataFile,
                                      const std::string& section )
{
    M_blockPrecs->setDataFromGetPot(  dataFile,
                                      section ,
                                      M_recompute.size());
}


void ComposedDN::coupler(map_shared_ptrtype& map,
                         const std::map<ID, ID>& locDofMap,
                         const vector_ptrtype& numerationInterface,
                         const Real& timeStep)
{
    UInt totalDofs( map->getMap(Unique)->NumGlobalElements() );
    UInt solidAndFluid(M_offset[solid]+1+M_FESpace[solid]->map().getMap(Unique)->NumGlobalElements());

    matrix_ptrtype coupling(new matrix_type(*map));
    couplingMatrix( coupling,  (*M_couplingFlags)[solid], M_FESpace, M_offset, locDofMap, numerationInterface, timeStep);
    coupling->insertValueDiagonal( 1., M_offset[fluid]+1, M_offset[solid]+1 );
    coupling->insertValueDiagonal( 1., solidAndFluid, totalDofs+1 );
    M_coupling.push_back(coupling);

    coupling.reset(new matrix_type(*map));
    couplingMatrix( coupling,  (*M_couplingFlags)[fluid], M_FESpace, M_offset, locDofMap, numerationInterface, timeStep);
    coupling->insertValueDiagonal( 1. , M_offset[solid], solidAndFluid );
    coupling->insertValueDiagonal( 1. , solidAndFluid + nDimensions*numerationInterface->getMap().getMap(Unique)->NumGlobalElements(), totalDofs +1 );
    M_coupling.push_back(coupling);

}

int ComposedDN::solveSystem( const vector_type& rhs, vector_type& step, solver_ptrtype& linearSolver )
{
    assert(M_blockPrecs.get());

    if(!set())
    {
        for(UInt k=0; k < M_blocks.size(); ++k)
            push_back_precs(M_blocks[(*M_blockReordering)[k]]);
    }
    else
    {
        for(UInt k=0; k < M_blocks.size(); ++k)
        {
            if(M_recompute[(*M_blockReordering)[k]])
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


void ComposedDN::push_back_precs( matrix_ptrtype& Mat )
{
    M_blockPrecs->push_back(Mat);
}

void ComposedDN::replace_precs(  matrix_ptrtype& Mat, UInt position )
{
    M_blockPrecs->replace(Mat, position);
}


} // Namespace LifeV
