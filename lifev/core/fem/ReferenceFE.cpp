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
    @brief Base structure for a reference finite element

    @author Jean-Frederic Gerbeau
    @date 00-04-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#include <lifev/core/fem/ReferenceFE.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

ReferenceFE::ReferenceFE ( std::string name, FE_TYPE type, ReferenceShapes shape,
                           Int nbDofPerVertex, Int nbDofPerEdge, Int nbDofPerFace,
                           Int nbDofPerVolume, Int nbDof, Int nbLocalCoor, Int FEDim, const function_Type* phi,
                           const function_Type* dPhi, const function_Type* d2Phi, const function_Type* divPhi , const Real* refCoor,
                           DofPatternType patternType,
                           const ReferenceFE* bdRefFE ) :
    ReferenceElement ( name, shape, nbDof, nbLocalCoor, FEDim, phi, dPhi, d2Phi, divPhi, refCoor ),
    DOFLocalPattern ( nbDof, nbDofPerVertex, nbDofPerEdge, nbDofPerFace, nbDofPerVolume, patternType, nbLocalCoor ),
    M_boundaryFE ( bdRefFE ), M_type ( type )
{}

ReferenceFE::~ReferenceFE()
{}



}
