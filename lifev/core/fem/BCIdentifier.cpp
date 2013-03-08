//@HEADER
/*
************************************************************************

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

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Implementations for BCIdentifier.hpp

    @date 01-07-2002
    @author M. A. Fernandez
    @contributor Luca Bertagna <lbertag@emory.edu>
 */

#include <lifev/core/fem/BCIdentifier.hpp>

namespace LifeV
{

////////////////////
// BCIdentifierBase //
////////////////////

// ============================================= //
//          Constructors & Destructor            //
// ============================================= //

BCIdentifierBase::BCIdentifierBase (const BCIdentifierBase& id) : M_id ( id.M_id )
{
    // Nothing to be done here
}

// ============================= //
//            Methods            //
// ============================= //

void BCIdentifierBase::showMe ( std::ostream& output) const
{
    output << "\n Node id: " << M_id << '\n';
}

/////////////////////////
// BCIdentifierEssential //
/////////////////////////

// ============================================= //
//          Constructors & Destructor            //
// ============================================= //

BCIdentifierEssential::BCIdentifierEssential ( BCIdentifierEssential const& id ) : BCIdentifierBase ( id ),
    M_x ( id.M_x ),
    M_y ( id.M_y ),
    M_z ( id.M_z )
{
    // Nothing to be done here
}

// ============================= //
//            Methods            //
// ============================= //

void BCIdentifierEssential::showMe ( std::ostream& output) const
{
    output << "\nNode id:" << M_id << '\n';
    output << "Node coordinates:\n";
    output << "              x = " << M_x;
    output << "              y = " << M_y;
    output << "              z = " << M_z;
}

///////////////////////
// BCIdentifierNatural //
///////////////////////

// ============================================= //
//          Constructors & Destructor            //
// ============================================= //

BCIdentifierNatural::BCIdentifierNatural ( const ID& id, const std::vector<ID>& localToGlobal ) : BCIdentifierBase ( id )
{
    M_localToGlobal.reserve ( localToGlobal.size() );
    M_localToGlobal.insert ( M_localToGlobal.end(), localToGlobal.begin(), localToGlobal.end() );
}

BCIdentifierNatural::BCIdentifierNatural ( const ID& id ) : BCIdentifierBase ( id )
{
    // Nothing to be done here
}

BCIdentifierNatural::BCIdentifierNatural ( BCIdentifierNatural const& id ) : BCIdentifierBase ( id ),
    M_localToGlobal ( id.M_localToGlobal )
{
    // Nothing to be done here
}

// ============================= //
//            Methods            //
// ============================= //

void BCIdentifierNatural::showMe ( std::ostream& output) const
{
    output << "\nNode id:" << M_id << '\n';
    output << "Local-to-global map:\n";

    int i (0);
    for ( std::vector<ID>::const_iterator it = M_localToGlobal.begin(); it != M_localToGlobal.end(); ++it, ++i )
    {
        output << "Local id: " << i << "  -->  Global id: " << *it << '\n';
    }
}

} // Namespace LifeV
