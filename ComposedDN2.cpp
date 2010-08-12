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
    @brief A short description of the file content

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 14 Jun 2010

    A more detailed description of the file (if necessary)
 */

#include <ComposedDN2.hpp>

namespace LifeV {

void ComposedDN2::coupler(map_shared_ptrtype map,
                          const std::map<ID, ID>& locDofMap,
                          const vector_ptrtype numerationInterface,
                          const Real& timeStep)
{
    UInt one(1.);

    matrix_ptrtype coupling2(new matrix_type(*map));
    //coupling.reset(new matrix_type(*map));
    couplingMatrix( coupling2, 8, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep );
    coupling2->insertValueDiagonal( one, 1 ,M_offset[0]+1);
    coupling2->insertValueDiagonal( one, M_offset[0]+1+M_FESpace[0]->map().getMap(Unique)->NumGlobalElements(), map->getMap(Unique)->NumGlobalElements()+1);
    M_coupling.push_back(coupling2);

    matrix_ptrtype coupling(new matrix_type(*map));
    couplingMatrix( coupling,  6, M_FESpace, M_offset, locDofMap, numerationInterface, timeStep);
    coupling->insertValueDiagonal( one, M_FESpace[0]->map() , M_offset[0]);/*dFESpace*/
    M_coupling.push_back(coupling);
    //blockView->insertValueDiagonal( one, *M_interfaceMap, offset2+M_dMap->getMap(Unique)->NumGlobalElements());
    //super::super::blockAssembling();
}

void ComposedDN2::blockAssembling()
{
    swap(0,1);
    super::blockAssembling();
}

void ComposedDN2::replace_matrix( const matrix_ptrtype& oper, UInt position)
{
    super::replace_matrix(oper, position);
}


void ComposedDN2::replace_precs( matrix_ptrtype& oper, UInt position)
{
    super::replace_precs(oper, position);
}


} // Namespace LifeV
