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
    @brief A short description of the file content

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 10 May 2010

    A more detailed description of the file (if necessary)
 */

#include <refFEScalar.hpp>

namespace LifeV {

RefFEScalar::RefFEScalar( std::string name, FE_TYPE type, ReferenceShapes shape,
                          int nbDofPerVertex, int nbDofPerEdge, int nbDofPerFace,
                          int nbDofPerVolume, int nbDof, int nbCoor, const Fct* phi,
                          const Fct* dPhi, const Fct* d2Phi, const Real* refCoor,
                          DofPatternType patternType,
                          const RefFE* bdRefFE, const ValuesToValuesFct nodalToFE ) :
    RefFE( name, type, shape,nbDofPerVertex,nbDofPerEdge,nbDofPerFace,
           nbDofPerVolume, nbDof, nbCoor,1, phi, dPhi, d2Phi, static_cast<Fct*>(NULL), refCoor,
           patternType, bdRefFE ),
    M_nodalToFEValues(nodalToFE)
{
    CONSTRUCTOR( "RefScalar" );
}

} // Namespace LifeV
