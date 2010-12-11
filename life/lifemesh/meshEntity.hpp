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
    @brief This file contains the MeshEntity and MeshEntityWithBoundary classes.

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @maintainer Tiziano Passerini <tiziano@mathcs.emory.edu>

    The classes included in this file are useful to store the identifiers of
    the different structures stored in the meshes.
 */

#ifndef MESHENTITY_H
#define MESHENTITY_H 1

#include <life/lifecore/life.hpp>

namespace LifeV
{
//! This is the base class to store the identifiers.
/*!
   In this class, there are two identifiers stored:
    <ol>
        <li> The global identifier is defined for the whole mesh.
        <li> The local identifier is defined as the identifier for
             a particular processor.
    </ol>

    When running the code in serial (1 processor), the identifiers
    are then the same.

    This class provides the method and operators to handle easily
    the identifiers.

    Note: Documentation by Samuel Quinodoz, implementation anterior
    to the documentation, without name of author.
 */
class MeshEntity
{
public:
    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    /*!
       Using this constructor, both identifiers are set to 0.
     */
    MeshEntity():
            M_id( 0 ),
            M_localId( 0 )
    {};

    //! Copy Constructor
    MeshEntity(const MeshEntity& meshEntity):
            M_id     (meshEntity.M_id),
            M_localId(meshEntity.M_localId)
    {};

    //! Constructor with a single value for both identifiers.
    /*!
       @param id The value for both identifers.
     */
    MeshEntity( ID id ):
            M_id( id ),
            M_localId( id )
    {};

    //! Full constructor, where both identifiers are specified.
    /*!
       @param id The value for the global ID.
       @param lid The value for the local ID.
     */
    MeshEntity( ID id, ID lid ):
            M_id( id ),
            M_localId( lid )
    {};

    //! Destructor
    ~MeshEntity()
    {};

    //@}


    //! @name Methods
    //@{


    //! Displays the informations stored by this class
    void showMe( std::ostream& output = std::cout ) const
    {
        output << " Global ID : " << M_id << " -- " << " Local ID " << M_localId << std::endl;
    };


    //@}


    //! @name Operators
    //@{

    //! Equivalence operator that checks if BOTH identifiers are the same.
    /*!
      @param e The mesh entity to be compared with.
     */
    bool operator==(const MeshEntity & e ) const
    {
        bool res = ( ( M_id == e.id() ) && ( M_localId == e.M_localId ));
        return res;
    };

    //! Relation operator that performs the same comparison on the GLOBAL identifier.
    /*!
      @param e The mesh entity to be compared with.
    */
    bool operator<=(const MeshEntity & e ) const
    {
        return M_id <= e.id();
    };

    //! Relation operator that performs the same comparison on the GLOBAL identifier.
    /*!
      @param e The mesh entity to be compared with.
    */
    bool operator>=(const MeshEntity & e ) const
    {
        return M_id >= e.id();
    };

    //@}


    //! @name Set Methods
    //@{

    //! Method to set the global identifier.
    /*!
      @param id The new global identifier.
    */
    inline void setId( ID id)
    {
        M_id = id;
    };

    //! Method to set the local identifier.
    /*!
      @param id The new local identifier.
    */
    inline void setLocalId( ID id)
    {
        M_localId = id;
    };

    //@}


    //! @name Get Methods
    //@{

    //! Method to get the global identifier.
    /*!
      @return The global identifier.
    */
    inline const ID & id() const
    {
        return M_id;
    };

    //! Method to get the local identifier.
    /*!
      @return The local identifier.
    */
    inline const ID & localId() const
    {
        return M_localId;
    };

    //@}

private:
    ID M_id;
    ID M_localId;
};



//! MeshEntityWithBoundary - This is a MeshEntity with an additional information on the boundary.
/*!
  The additional boolean that is stored is used to know if the MeshEntity is on the boundary
  or not. This class provides all the methods needed to handle easily this additional information.

  Note: Documentation by Samuel Quinodoz, implementation anterior
  to the documentation, without name of author.
 */
class MeshEntityWithBoundary : public MeshEntity
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    /*!
      This constructor calls the empty constructor of MeshEntity and
      sets the boundary indicator to false.
    */
    MeshEntityWithBoundary() : MeshEntity(), M_boundary( false )
    {};

    //! Copy constructor
    MeshEntityWithBoundary( const MeshEntityWithBoundary& meshEntityWithBoundary ) :
            MeshEntity( meshEntityWithBoundary ),
            M_boundary( meshEntityWithBoundary.M_boundary )
    {};

    //! Specific constructor
    /*!
      This is the "full" constructor for this class.
      @param id The identifier to be set for both (global and local) identifiers. Use
      a set method if you want different identifiers.
      @param boundary The value of the boundary indicator.
    */
    MeshEntityWithBoundary( ID id, bool boundary = false ) :
            MeshEntity( id ),
            M_boundary( boundary )
    {};

    //! Destructor
    ~MeshEntityWithBoundary()
    {};

    //@}


    //! @name Methods
    //@{


    //! Display the informations stored by this class
    void showMe( std::ostream& output = std::cout ) const
    {
        output << " Global ID : " << id() << " -- " << " Local ID " << localId();
        if (M_boundary)
        {
            output << " -- Boundary ";
        };
        output << std::endl;
    };


    //@}


    //! @name Set Methods
    //@{


    //! Set method for the boundary indicator
    /*!
      @param boundary The value to be set for the boundary indicator.
    */
    inline void setBoundary(const bool& boundary)
    {
        M_boundary = boundary;
    };


    //@}

    //! @name Get Methods
    //@{


    //! Tells if it is on the boundary
    inline const bool & boundary() const
    {
        return M_boundary;
    };


    //@}

private:
    bool M_boundary;
};


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
#endif
