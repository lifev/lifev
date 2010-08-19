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
    @date 09 Aug 2010
 */

#include <ComposedDND.hpp>

namespace LifeV {

void ComposedDND::blockAssembling( )
{
    if(!M_swapped)
    {
        super::super::super::swap(M_coupling[2], M_coupling[3]);
        M_swapped = true;
    }

    addToCoupling(M_coupling[3], 2);

    super::blockAssembling();
}

} // Namespace LifeV
