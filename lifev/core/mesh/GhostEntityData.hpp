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
    @file
    @brief Ghost entity data structure

    @author Radu Popescu <radu.popescu@epfl.ch>
    @date 17-11-2011
 */

#ifndef GHOST_ENTITY_DATA_H
#define GHOST_ENTITY_DATA_H 1

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

//! Ghost entity data structure
/*!
    @author Radu Popescu <radu.popescu@epfl.ch>

    TODO: add description
 */
struct GhostEntityData
{
    //! ID in the current sub-domain of the face.
    ID localFacetId;

    //! ID of element that faces the local one on the other sub-domain.
    ID ghostElementLocalId;

    //! Position on the ghost element.
    UInt ghostElementPosition;

    friend std::ostream& operator<< (std::ostream& out, GhostEntityData const& ged);
};

} // Namespace LifeV

#endif // GHOST_ENTITY_DATA_H
