/*
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
    MeshEntity():
        M_id( 0 ),
        M_localId( 0 )
        {}
    ;
    MeshEntity(const MeshEntity& meshEntity):
        M_id     (meshEntity.M_id),
        M_localId(meshEntity.M_localId)
        {}
    ;
    MeshEntity( ID id ):
        M_id( id ),
        M_localId( id )
        {}

    MeshEntity( ID id, ID lid ):
        M_id( id ),
        M_localId( lid )
        {}
    ;
    ID id() const
        {
            return M_id;
        }

    ID localId() const
        {
            return M_localId;
        }

    void setId( ID id)
        {
            M_id = id;
        }

    void setLocalId( ID id)
        {
            M_localId = id;
        }



//     ID & id()
//     {
//         return M_id;
//     }

    bool operator==( MeshEntity & e ) const
    {
        bool res = ( ( M_id == e.id() ) && ( M_localId == e.M_localId ));
        return res;
    };


    bool operator<=( MeshEntity & e ) const
    {
        return M_id <= e.id();
    };

    bool operator>=( MeshEntity & e ) const
    {
        return M_id >= e.id();
    };

private:
    ID M_id;
    ID M_localId;

};


//! Base class with boundary
/*! Contains info on boundary position */
class MeshEntityWithBoundary : public MeshEntity
{
public:
    MeshEntityWithBoundary() : MeshEntity(), _boundary( false )
    {}
    ;
    MeshEntityWithBoundary( const MeshEntityWithBoundary& meshEntityWithBoundary ) :
        MeshEntity( meshEntityWithBoundary ),
        _boundary( meshEntityWithBoundary._boundary )
    {}
    ;
    MeshEntityWithBoundary( ID i, bool boundary = false ) :
        MeshEntity( i ),
        _boundary( boundary )
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
    OrientedMeshEntity() :
        MeshEntity(),
        _orient( true )
    {}
    ;
    OrientedMeshEntity(const OrientedMeshEntity& orientedMeshEntity) :
        MeshEntity( orientedMeshEntity ),
        _orient(orientedMeshEntity._orient )
    {}
    ;
    OrientedMeshEntity( ID i, bool o = true ) :
        MeshEntity( i ),
        _orient( o )
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
