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
 * @brief Implementation of the file refEle.hpp
 */

#include <life/lifefem/refEle.hpp>

namespace LifeV
{

RefEle::RefEle( std::string name, ReferenceShapes shape, UInt nbDof, UInt nbCoor, UInt FEDim,
                const Fct* phi, const Fct* dPhi, const Fct* d2Phi, const Fct* divPhi, const Real* refCoor ) :
    M_phi( phi ),
    M_dPhi( dPhi ),
    M_d2Phi( d2Phi ),
    M_divPhi( divPhi),
    M_refCoor( refCoor ), 
    
    M_name( name ),
    M_shape( shape ),
    M_nbDof( nbDof ),
    M_nbCoor( nbCoor ),
    M_FEDim( FEDim )
{
    CONSTRUCTOR( "RefEle" );
}

RefEle::~RefEle()
{
    DESTRUCTOR( "RefEle" );
}
    
std::vector<GeoVector> 
RefEle::refCoor() const
{
    std::vector<GeoVector> coordinates(M_nbDof, GeoVector(3,0));
    for (UInt i(0); i<M_nbDof; ++i)
    {
        coordinates[i][0]=M_refCoor[3*i];
        coordinates[i][1]=M_refCoor[3*i+1];
        coordinates[i][2]=M_refCoor[3*i+2];
    }
    return coordinates;
}


} // Namespace LifeV
