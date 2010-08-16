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
    @date 14 Jun 2010
 */

#include <ComposedDN2.hpp>

namespace LifeV {

void ComposedDN2::coupler(map_shared_ptrtype& map,
                          const std::map<ID, ID>& locDofMap,
                          const vector_ptrtype& numerationInterface,
                          const Real& timeStep)
{
    super::coupler( map, locDofMap, numerationInterface, timeStep );
    M_blockReordering.resize(3);
    M_blockReordering[0] = fluid;
    M_blockReordering[1] = solid;
    M_blockReordering[2] = mesh;

}

} // Namespace LifeV
