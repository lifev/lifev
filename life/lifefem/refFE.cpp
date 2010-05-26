/*-*- mode: c++ -*-
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <life/lifefem/refFE.hpp>

namespace LifeV
{

RefFE::RefFE( std::string name, FE_TYPE type, ReferenceShapes shape,
              int nbDofPerVertex, int nbDofPerEdge, int nbDofPerFace,
                      int nbDofPerVolume, int nbDof, int nbCoor, int FEDim, const Fct* phi,
              const Fct* dPhi, const Fct* d2Phi, const Fct* divPhi , const Real* refCoor,
              DofPatternType patternType,
              const RefFE* bdRefFE ) :
    RefEle( name, shape, nbDof, nbCoor, FEDim, phi, dPhi, d2Phi, divPhi, refCoor ),
    LocalDofPattern( nbDof, nbDofPerVertex, nbDofPerEdge, nbDofPerFace, nbDofPerVolume, patternType ),
    M_boundaryFE( bdRefFE ), M_type( type )
{
    CONSTRUCTOR( "RefFE" );
}

RefFE::~RefFE()
{
    DESTRUCTOR( "RefFE" );
}



}
