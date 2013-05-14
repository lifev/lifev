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
    @brief This file contains the MeshEntity class.

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @contributor Antonio Cervone <ant.cervone@gmail.com>
    @maintainer Tiziano Passerini <tiziano@mathcs.emory.edu>

    The classes included in this file are useful to store the identifiers of
    the different structures stored in the meshes.
 */

#ifndef MESHENTITY_H
#define MESHENTITY_H 1

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

//! available bit-flags for different geometric properties
namespace EntityFlags
{
const flag_Type DEFAULT             ( 0x00 );
const flag_Type PHYSICAL_BOUNDARY   ( 0x01 );
const flag_Type INTERNAL_INTERFACE  ( 0x02 );
const flag_Type SUBDOMAIN_INTERFACE ( 0x04 );
const flag_Type OVERLAP             ( 0x08 );
const flag_Type CUTTED              ( 0x10 );
const flag_Type VERTEX              ( 0x20 );
const flag_Type GHOST               ( 0x40 );
// @note remember to update ALL value in order to encompass all flags
const flag_Type ALL                 ( 0x7F );

const UInt number                   (    7 );

std::string name ( const flag_Type& flag );

}// namespace EntityFlags

//! This is the base class to store basic properties of any mesh entity
/*!
 * The basic properties of a mesh entity is to have identifiers and flags.
 * Identifiers give information about the numbering of the entity (typically in a mesh)
 * flags are attributes af the entity, like being on the boundary.
 *
   In this class, there are two identifiers stored:
    <ol>
        <li> The global identifier, normally used to uniquely identify the
         entity in the whole mesh.
        <li> The local identifier, which is usually unique to the entity on a given
        partition. If the entity is stored into a
        meshEntityContainer, then localId must be the position
        of the entity in the container.
    </ol>

    When running the code in serial (1 processor), the identifiers
    are then the same (but this is not necessary)

    This class provides the method and operators to handle local and global
    identifiers easily.

    Note: Documentation by Samuel Quinodoz, implementation anterior
    to the documentation: Luca Formaggia

  The additional fag that is stored is used to specify geometrical properties of the entity
  such as the ones specified above in the EntityFlags namespace.
  @todo The marker ID should be taken away form marker class and added to the MeshEntity directly
 */

class MeshEntity
{
public:
    //! Indicator for local or global id
    /*!
     * This enum helps to develop methods to operate on local or globalID
     *
     */
    enum SwitchId {LOCALID = 0, GLOBALID = 1};
    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    /*!
       Using this constructor, both identifiers are set to NotAnId.
     */
    MeshEntity() :
        M_id ( NotAnId ),
        M_localId ( NotAnId ),
        M_flag ( EntityFlags::DEFAULT )
    {}

    //! Constructor with a single value for both identifiers.
    /*!
       @param id The value for both identifers.
       @param flag The value of the flag to assign (optional)
     */
    MeshEntity ( const ID& id, const flag_Type& flag = EntityFlags::DEFAULT ) :
        M_id ( id ),
        M_localId ( id ),
        M_flag ( flag )
    {}

    //! Full constructor, where both identifiers are specified.
    /*!
       @param id The value for the global ID.
       @param lid The value for the local ID.
       @param flag The value of the flag to assign (optional)
     */
    MeshEntity ( const ID& id, const ID& lid, const flag_Type& flag = EntityFlags::DEFAULT ) :
        M_id ( id ),
        M_localId ( lid ),
        M_flag ( flag )
    {}

    //! backward-compatible constructor
    /*!
      @param id The identifier to be set for both (global and local) identifiers. Use
      a set method if you want different identifiers.
      @param boundary The value of the boundary indicator.
    */
    MeshEntity ( const ID& id, const bool& boundary ) :
        M_id ( id ),
        M_localId ( id )
    {
        // NOTE: this is conservative, if PHYSICAL_BOUNDARY = 0x01
        // this can be done as
        // M_flag = static_cast<flag_Type> boundary;
        M_flag = boundary ? EntityFlags::PHYSICAL_BOUNDARY : EntityFlags::DEFAULT;
    }

    //! Destructor
    virtual ~MeshEntity() {}

    //@}

    //! @name Methods
    //@{
    //! Displays the informations stored by this class
    void showMe ( std::ostream& output = std::cout ) const;
    //@}


    //! @name Set Methods
    //@{

    //! Method to set the global identifier.
    /*!
     * @todo change the name in setGlobalId
      @param id The new global identifier.
    */
    inline void setId ( const ID& id)
    {
        M_id = id;
    }

    //! Method to set the local identifier.
    /*!
      @param id The new local identifier.
    */
    inline void setLocalId ( const ID& id)
    {
        M_localId = id;
    }

    //! Set method for the boundary indicator
    /*!
      @param boundary The value to be set for the boundary indicator.
    */
    void setBoundary (const bool& boundary)
    {
        if ( boundary )
        {
            M_flag = Flag::turnOn  ( M_flag, EntityFlags::PHYSICAL_BOUNDARY );
        }
        else
        {
            M_flag = Flag::turnOff ( M_flag, EntityFlags::PHYSICAL_BOUNDARY );
        }
    }

    //! Replace method for the entity flag
    /*!
      @param flag The value to be set for the entity flag.
      @note Beware, it sets the entire flag; if you want to add a flag use | or setFlag
    */
    void replaceFlag ( const flag_Type& flag )
    {
        M_flag = flag;
    }

    //! Sets a flag
    /**
     * @param flag The flag to be added
     */
    void setFlag ( const flag_Type& flag )
    {
        M_flag = Flag::turnOn (flag, M_flag);
    }
    //! Remove a flag
    /**
     * @param flag The flag to be removed
     */
    void unSetFlag ( const flag_Type& flag )
    {
        M_flag = Flag::turnOff (M_flag, flag);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Method to get the global identifier.
    /*!
      @return The global identifier.
    */
    inline const ID& id() const
    {
        return M_id;
    }

    //! Method to get the local identifier.
    /*!
      @return The local identifier.
    */
    inline const ID& localId() const
    {
        return M_localId;
    }

    //! Tells if it is on the boundary
    bool boundary() const
    {
        return Flag::testOneSet ( M_flag, EntityFlags::PHYSICAL_BOUNDARY );
    }

    //! Tells if the entity is owned by current process
    bool isOwned() const
    {
        return Flag::testOneNotSet ( M_flag, EntityFlags::GHOST );
    }

    //! returns the entity flag
    const flag_Type& flag() const
    {
        return M_flag;
    }

    //@}

private:
    ID M_id;
    ID M_localId;
    flag_Type M_flag;
};

namespace MeshEntityUtility
{
/*! @defgroup MeshEntityUtilities
 *  Utilities to get local or global ID according to a switch
 *  The template parameter is in fact a MeshEntity::SwitchId
 *  enumerator which takes values MeshEntity::LOCALID or
 *  MeshEntity::GLOBALID.
 *  Getters are also implemented as functors for efficiency reason
 *  (allow inlining)
 *  @todo Go to a separate file
 * @{
 */
template<int Selector>
inline ID getID (MeshEntity const&);
//! Generic definition of setter
template<int Selector>
inline void setID (MeshEntity&, const ID);

//! Specialization for global id
template<>
inline ID getID<MeshEntity::GLOBALID> (MeshEntity const& entity)
{
    return entity.id();
}

//! Specialization for local id
template<>
inline ID getID<MeshEntity::LOCALID> (MeshEntity const& entity)
{
    return entity.localId();
}

//! Specialization for global id
template<>
inline void setID<MeshEntity::GLOBALID> (MeshEntity& entity, const ID id)
{
    entity.setId (id);
}

//! Specialization for local id
template<>
inline void setID<MeshEntity::LOCALID> (MeshEntity& entity, const ID id)
{
    entity.setLocalId (id);
}
//! Generic definition of the functor to extract the local or global ID
/*
 * The functor is useful for allowing inlining and use it with
 * std compliant containers.
 */
template<int Selector>
struct IdGetter
{
    inline ID operator() (MeshEntity const& entity) const
    {
        return getID<Selector> (entity);
    };
};



/** @} */ // end of MeshEntityUtilities
} // end of namespace MeshEntityUtility

//! Mesh Entity with an orientation
/*! \note Not Used so far! */
/*class OrientedMeshEntity: public MeshEntity
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
    }; */


} // End of namespace LifeV

#endif /* MESHENTITY_H */
