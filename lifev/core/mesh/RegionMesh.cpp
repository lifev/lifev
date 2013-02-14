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

/*!
 *  @file
 *  @brief File containing 3D Mesh Class Implementation
 *
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *  @author Miguel Fernandez
 *
 *  @contributor Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 *  @mantainer Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 */

#include <lifev/core/util/Switch.hpp>

namespace LifeV
{
void set_switches_for_regionmesh ( Switch& sw )
{
    sw.create ( "HAS_ALL_FACETS" );
    sw.create ( "HAS_ALL_RIDGES" );
    sw.create ( "HAS_BOUNDARY_FACETS" );
    sw.create ( "HAS_BOUNDARY_RIDGES" );
    sw.create ( "HAS_ELEMENT_TO_FACETS" );
    sw.create ( "HAS_ELEMENT_TO_RIDGES" );
    sw.create ( "HAS_BEEN_CHECKED" );
    sw.create ( "FACETS_HAVE_ADIACENCY" );
}
}
