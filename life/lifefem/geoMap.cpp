//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

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
 * @file
 * @brief Implementation of the file geoMap.hpp
 */

#include <life/lifefem/geoMap.hpp>

namespace LifeV
{

GeoMap::GeoMap( std::string _name, ReferenceShapes _shape,
                UInt _nbDof, UInt _nbCoor,
                const Fct* phi, const Fct* dPhi, const Fct* d2Phi,
                const Real* _refCoor,
                const GeoMap* bdMap ) :
    RefEle( _name, _shape, _nbDof, _nbCoor,1, phi, dPhi, d2Phi, static_cast<Fct*>(NULL),  _refCoor ),
        M_boundaryMap( bdMap )
{
    CONSTRUCTOR( "GeoMap" );
}
GeoMap::~GeoMap()
{
    DESTRUCTOR( "GeoMap" );
}

}
