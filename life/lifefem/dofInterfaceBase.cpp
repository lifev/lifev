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
    @brief Base Class for interfacing dofs between two meshes

    @author M.A. Fernandez and V. Martin
    @date 00-11-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    This file contains the class which may be used to update and hold the connections between the dof
    on two matching meshes.
 */

#include <life/lifefem/dofInterfaceBase.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

DofInterfaceBase::DofInterfaceBase()
{}

// ===================================================
// Methods
// ===================================================


ID DofInterfaceBase::getInterfaceDof( const ID& i ) const
{
    std::map<ID, ID>::const_iterator it = M_locDofMap.find( i );
    if ( it == M_locDofMap.end() )
    {
        std::cout << i << " : " << std::flush;
        ERROR_MSG( "DofInterfaceBase::getInterfaceDof : Dof number not found" );
    }
    return it->second;
}

bool DofInterfaceBase::isMyInterfaceDof( const ID& i ) const
{
    std::map<ID, ID>::const_iterator it = M_locDofMap.find( i );
    return ( it != M_locDofMap.end() );
}

std::ostream& DofInterfaceBase::showMe( bool verbose, std::ostream& out ) const
{
    out << "------------------------------" << std::endl;
    out << "\tNumber of Dof connections (M_locDofMap):" << M_locDofMap.size() << std::endl;
    if ( verbose )
    {
        UInt count( 0 ), lines( 10 );
        out << "List of connections between Dof: (global, local)";
        for ( std::map<ID, ID>::const_iterator it = M_locDofMap.begin(); it != M_locDofMap.end(); ++it )
        {
            if ( count++ % lines == 0 )
            {
                out << std::endl;
            }
            out << "(" << it->first << "," << it->second << ")\t";
        }
        out << std::endl;
    }
    out << "------------------------------" << std::endl;
    return out;
}

void DofInterfaceBase::buildInverse( const DofInterfaceBase& dofBase)
{
    for ( std::map<ID, ID>::const_iterator it = dofBase.M_locDofMap.begin(); it != dofBase.M_locDofMap.end(); ++it )
    {
        M_locDofMap[it->second] = it->first;
    }

}

// ===================================================
// Set Methods
// ===================================================

ID DofInterfaceBase::nbInterfaceDof() const
{
    return M_locDofMap.size();
}

// ===================================================
// Get Methods
// ===================================================

void DofInterfaceBase::set( ID key, ID value )
{
    M_locDofMap.find(key)->second = value;
}

}
