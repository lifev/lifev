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
    @brief Base class for RefFE and GeoMap

    @author Jean-Frederic Gerbeau
    @date 00-04-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#include <life/lifefem/refEle.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

RefEle::RefEle( std::string name, ReferenceShapes shape, UInt nbDof, UInt nbCoor, UInt feDim,
                const function_Type* phi, const function_Type* dPhi, const function_Type* d2Phi,
                const function_Type* divPhi, const Real* refCoor ) :
        M_phi( phi ),
        M_dPhi( dPhi ),
        M_d2Phi( d2Phi ),
        M_divPhi( divPhi),
        M_refCoor( refCoor ),

        M_name( name ),
        M_shape( shape ),
        M_nbDof( nbDof ),
        M_nbCoor( nbCoor ),
        M_feDim( feDim )
{
}

RefEle::~RefEle()
{
}

// ===================================================
// Methods
// ===================================================

std::vector<GeoVector>
RefEle::refCoor() const
{
    std::vector<GeoVector> coordinates(M_nbDof, GeoVector(3));
    for (UInt i(0); i<M_nbDof; ++i)
    {
        coordinates[i][0]=M_refCoor[3*i];
        coordinates[i][1]=M_refCoor[3*i+1];
        coordinates[i][2]=M_refCoor[3*i+2];
    }
    return coordinates;
}


} // Namespace LifeV
