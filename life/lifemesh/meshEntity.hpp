/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
/*! file MeshEntity.h */
#ifndef _MESHENTITY_HH_
#define _MESHENTITY_HH_
#include <life/lifecore/life.hpp>

namespace LifeV
{
//using namespace std; //Pretty useless

//! Base class of all Mesh Entities
/*! It contains the Entity ID and it defined the comparison operators */
class MeshEntity
{
public:
    MeshEntity() : _id( 0 )
    {}
    ;
    MeshEntity( ID i ) : _id( i )
    {}
    ;
    ID id() const
    {
        return _id;
    }
    ID & id()
    {
        return _id;
    }
    bool operator==( MeshEntity & e ) const
    {
        return _id == e.id();
    };
    bool operator<=( MeshEntity & e ) const
    {
        return _id <= e.id();
    };
    bool operator>=( MeshEntity & e ) const
    {
        return _id >= e.id();
    };
protected:
    ID _id;
};


//! Base class with boundary
/*! Contains info on boundary position */
class MeshEntityWithBoundary : public MeshEntity
{
public:
    MeshEntityWithBoundary() : MeshEntity(), _boundary( false )
    {}
    ;
    MeshEntityWithBoundary( ID i, bool boundary = false ) : MeshEntity( i ), _boundary( boundary )
    {}
    ;
    //! Tells if  item is on the boundary
    bool boundary() const
    {
        return _boundary;
    };
    //! Changes boundary indicator
    bool & boundary()
    {
        return _boundary;
    };
protected:
    bool _boundary;
};


//! Mesh Entity with an orientation
/*! \note Not Used so far! */
class OrientedMeshEntity: public MeshEntity
{
public:
    OrientedMeshEntity() : MeshEntity(), _orient( true )
    {}
    ;
    OrientedMeshEntity( ID i, bool o = true ) : MeshEntity( i ), _orient( o )
    {}
    //! Return entity orientation
    bool orientation() const
    {
        return _orient;
    }
    //! Assigns entity orientation
    bool & orientation()
    {
        return _orient;
    }
protected:
    bool _orient;
};
}
#endif
