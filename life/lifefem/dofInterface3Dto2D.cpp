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
    @brief Class for connecting the dof of a mesh (3D) and an interface (2D)
         that lives on the boundary of the mesh.

    @author Vincent Martin
    @date 00-02-2003

    This file contains the class which may be used to update and hold
    the connections between the dof of a mesh (3D) and an interface (2D)
    that lives on the boundary of the mesh. The interface is referenced
    by a flag.

    @contributor M.A. Fernandez
                 Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#include <life/lifefem/dofInterface3Dto2D.hpp>

namespace LifeV
{

// ===================================================
// Helper
// ===================================================


void RemoveMultiple( const std::list<ID> & listToTreat, std::list< std::pair<ID, ID> > & finalList )
{
    ID counter = 1;
    std::list<ID> temporaryList( listToTreat );

    //! Sort the list
    temporaryList.sort();

    //! initialize the new list
    std::pair <ID, ID> p0( temporaryList.front() , counter );
    finalList.push_back( p0 );

    //! We remove the multiple occurences :
    for ( std::list<ID>::iterator it = temporaryList.begin() ; it != temporaryList.end() ; ++ it )
    {
        if ( ( *it ) != finalList.back().first )
        {
            counter ++ ;
            //! Add to the list the new value
            std::pair <ID, ID> p( ( *it ) , counter );
            finalList.push_back( p );
        }
    }
}

// ===================================================
// Constructors & Destructor
// ===================================================

DofInterface3Dto2D::DofInterface3Dto2D( const LocalDofPattern& refFE, const Dof& dof1 ) :
        M_interfaceFlag( 0 ), M_refFE1( &refFE ), M_dof1( &dof1 )
{
    M_finalized = false;
}

// ===================================================
// Methods
// ===================================================

void
DofInterface3Dto2D::setup( const LocalDofPattern& refFE1, const Dof& dof1 )
{
    M_refFE1    = &refFE1;
    M_dof1      = &dof1;

    M_finalized = false;
}

void DofInterface3Dto2D::clearLists()
{
    M_vertexPerFaceList.clear();
    M_vertexList.clear();
}

std::ostream& DofInterface3Dto2D::showMe( bool verbose, std::ostream& out ) const
{
    out << "------------------------------" << std::endl;
    out << "myDofInterface reference: " << M_interfaceFlag << std::endl;
    out << "Number of face connections (M_faceList): " << M_faceList.size() << std::endl;
    if ( verbose )
    {
        unsigned int count( 0 ), lines( 10 );
        out << "\tList of connections between Faces: (global, local)";
        for ( std::vector< std::pair<ID, ID> >::const_iterator i = M_faceList.begin(); i != M_faceList.end(); ++i )
        {
            if ( count++ % lines == 0 )
            {
                out << std::endl;
            }
            out << "(" << i->first << "," << i->second << ")\t";
        }
        out << std::endl;
    }
    out << "Number of connections between Vertices (M_vertexList): " << M_vertexList.size() << std::endl;
    if ( verbose )
    {
        unsigned int count( 0 ), lines( 10 );
        out << "\tList of connections between Vertices: (global, local)";
        for ( std::list< std::pair<ID, ID> >::const_iterator it = M_vertexList.begin(); it != M_vertexList.end(); ++it )
        {
            if ( count++ % lines == 0 )
            {
                out << std::endl;
            }
            out << "(" << it->first << "," << it->second << ")\t";
        }
        out << std::endl;
    }
    //! print M_locDofMap
    showMe( verbose, out );

    out << "------------------------------" << std::endl;
    return out;
}

// ===================================================
// Operators
// ===================================================


ID DofInterface3Dto2D::operator[] ( const UInt& i ) const
{
    ASSERT_PRE( M_finalized, "The face List should be finalised before being accessed" );
    ASSERT_BD( i < M_faceList.size() );
    return M_faceList[ i ].first;  // M_faceList must be a vector!
}


DofInterface3Dto2D & DofInterface3Dto2D::operator=( const DofInterface3Dto2D& dofi )
{
    M_interfaceFlag = dofi.M_interfaceFlag;
    M_refFE1 = dofi.M_refFE1;
    M_dof1 = dofi.M_dof1;
    M_faceList = dofi.M_faceList;
    M_vertexPerFaceList = dofi.M_vertexPerFaceList; // (empty)
    M_vertexList = dofi.M_vertexList;
    M_localDofMap = dofi.M_localDofMap;
    M_finalized = dofi.M_finalized;

    return *this;
}


// ===================================================
// Private Methods
// ===================================================

ID DofInterface3Dto2D::vertex3Dto2D( const ID& idpoint3D ) const
{
    ASSERT_PRE( M_finalized, "The list of vertices must be finalized before accessing to the interface vertices." );
    for ( std::list< std::pair<ID, ID> >::const_iterator it = M_vertexList.begin(); it != M_vertexList.end(); ++it )
    {
        if ( it->first == idpoint3D )
            return it->second;
    }
    ERROR_MSG( "There is no such 3D index of vertex in the M_vertexList." );

    return ID();
}




}
