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

#include <refFEHybrid.hpp>

namespace LifeV {


// Costructor.
RefFEHybrid::RefFEHybrid( std::string name, FE_TYPE type, ReferenceShapes shape, UInt nbDofPerVertex, UInt nbDofPerEdge, 
                          UInt nbDofPerFace, UInt nbDofPerVolume, UInt nbDof, UInt nbCoor, const UInt& numberBoundaryFE, 
                          const StaticBdFE* boundaryFEList, const Real* refCoor, DofPatternType patternType ) :
        RefFE( name, type, shape, nbDofPerVertex, nbDofPerEdge, nbDofPerFace, nbDofPerVolume, 
                   nbDof, nbCoor,1,static_cast<Fct*>(NULL),  static_cast<Fct*>(NULL),
                   static_cast<Fct*>(NULL),  static_cast<Fct*>(NULL),refCoor, patternType, static_cast<RefFE*>(NULL)),
        M_numberBoundaryFE( numberBoundaryFE ), M_boundaryFEList( boundaryFEList )
{
    CONSTRUCTOR( "RefFEHybrid" );
}

// Destructor.
RefFEHybrid::~RefFEHybrid()
{
    DESTRUCTOR( "RefFEHybrid" );
}



} // Namespace LifeV
