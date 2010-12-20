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
    @brief CurrentFE living on the sides of the elements

    @author Jean-Frederic Gerbeau
    @date 00-09-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

 */

#include <life/lifefem/CurrentBoundaryFE.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

CurrentBdFE::CurrentBdFE( const RefFE& refFE, const GeoMap& geoMap, const QuadRule& qr ) :
        StaticBdFE( refFE, geoMap, qr )
{}

CurrentBdFE::CurrentBdFE( const RefFE& refFE, const GeoMap& geoMap ) :
        StaticBdFE( refFE, geoMap )
{}

CurrentBdFE::~CurrentBdFE()
{}

}
