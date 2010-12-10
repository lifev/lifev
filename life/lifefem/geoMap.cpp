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
    @brief This file contains the implementation of the GeoMap class (and an helper function)

    @author Jean-Frederic Gerbeau
    @date 00-04-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    This class contains the geometrical transformation that maps the reference
    element on the current element.
 */


#include <life/lifefem/geoMap.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

GeoMap::GeoMap( std::string name, ReferenceShapes shape,
                UInt nbDof, UInt nbCoor,
                const Fct* phi, const Fct* dPhi, const Fct* d2Phi,
                const Real* refCoor,
                const GeoMap* bdMap ) :
        RefEle( name, shape, nbDof, nbCoor,1, phi, dPhi, d2Phi, static_cast<Fct*>(NULL),  refCoor ),
        M_boundaryMap( bdMap )
{}
GeoMap::~GeoMap()
{}

}
