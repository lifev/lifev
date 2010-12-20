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

#include <MonolithicBlockComposedDND.hpp>

namespace LifeV
{

void MonolithicBlockComposedDND::blockAssembling( )
{
    if (!M_swapped)
    {
        super_Type::super_Type::super_Type::swap(M_coupling[2], M_coupling[3]);
        M_swapped = true;
    }

    addToCoupling(M_coupling[3], 2);

    super_Type::blockAssembling();
}

} // Namespace LifeV
