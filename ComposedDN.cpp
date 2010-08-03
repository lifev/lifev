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
                                      section );
}


void ComposedDN::coupler(map_shared_ptrtype map,
                         const std::map<ID, ID>& locDofMap,
                         const vector_ptrtype numerationInterface,
                         const Real& timeStep)
{
    UInt one(1.);
    matrix_ptrtype coupling(new matrix_type(*map));
    coupling->insertValueDiagonal( one, M_offset[1]+1, M_offset[0]+1);
    coupling->insertValueDiagonal( one, M_offset[0]+1+M_FESpace[0]->map().getMap(Unique)->NumGlobalElements(),  map->getMap(Unique)->NumGlobalElements()+1);
    coupling->GlobalAssemble();
    M_coupling.push_back(coupling);

    coupling.reset(new matrix_type(*map));
    couplingMatrix( coupling,  M_couplingFlag, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep) ;
    coupling->insertValueDiagonal( one, M_FESpace[0]->map()/*dFESpace*/ , M_offset[0] );
    coupling->GlobalAssemble();
    M_coupling.push_back(coupling);
}

int ComposedDN::solveSystem( const vector_type& rhs, vector_type& step, solver_ptrtype& linearSolver )
{
    assert(M_blockPrecs.get());

    if(!M_blockPrecs->getNumber())
    {
        for(UInt k=0; k < M_blocks.size(); ++k)
            push_back_precs(M_blocks[k]);
    }
    else
    {
        for(UInt k=0; k < M_blocks.size(); ++k)
        {
            if(M_recompute[k])
                replace_precs(M_blocks[k], k);
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
