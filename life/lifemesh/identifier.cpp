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
    @brief Implementations for identifier.hpp

    @date 01-07-2002
    @author M. A. Fernandez
    @contributor Luca Bertagna <lbertag@emory.edu>
 */

#include <life/lifemesh/identifier.hpp>

namespace LifeV
{

////////////////////
// IdentifierBase //
////////////////////

// ===================================================
// Constructors & Destructor
// ===================================================

IdentifierBase::IdentifierBase(const IdentifierBase & id)
{
    M_id = id.M_id;
}

/////////////////////////
// IdentifierEssential //
/////////////////////////

// ===================================================
// Constructors & Destructor
// ===================================================

IdentifierEssential::IdentifierEssential( IdentifierEssential const & id )
{
    M_id = id.M_id;
    M_x  = id.M_x;
    M_y  = id.M_y;
    M_z  = id.M_z;
}


///////////////////////
// IdentifierNatural //
///////////////////////

// ===================================================
// Constructors & Destructor
// ===================================================
IdentifierNatural::IdentifierNatural( const ID& id, const SimpleVect<ID>& localToGlobal ) : IdentifierBase( id )
{
    M_localToGlobal.reserve( localToGlobal.size() );
    M_localToGlobal.insert( M_localToGlobal.end(), localToGlobal.begin(), localToGlobal.end() );
}

IdentifierNatural::IdentifierNatural( const ID& id ) : IdentifierBase( id )
{
    // Nothing to be done here
}

IdentifierNatural::IdentifierNatural( IdentifierNatural const & id )
{
    M_id            = id.M_id;
    M_localToGlobal = id.M_localToGlobal;
}


} // Namespace LifeV
