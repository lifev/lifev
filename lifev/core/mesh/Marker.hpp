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
    @brief Basic definition of markers

    @contributor Luca Bertagna <lbertag@emory.edu>
    @date 00-00-0000

    Here we define the basic markers. Markers have two purposes:

    <ul>
        <li> To add an indicator (markerID_Type) to all geometry entities. In the
        base MarkerTrait this indicator is a long integer (aliased to
        EntityFLAG). The main purpose of the marker ID is to associate
        boundary conditions or material properties to the Geometry
        entity. </li>

        <li> There is a mechanism, via the MarkerCommon class template, to
        extend the mesh entities definitions by adding user defined data
        and methods. It this case the user should create its own marker
        class (maybe by inheriting from the Marker class, which
        provides the same interfaces of Marker + the user defined
        stuff.</li>
    </ul>

    The marker ID in the base class+some utilities to select between two marker IDs
    is provided by traits. In particular,

    <ul>
        <li> The MarkerIDStandardPolicy define the basic (compulsory) interface of
        any user defined Marker class.</li>

        <li> The Marker is a class template whose template argument is a
        MarkerTrait, defaulted to MarkerIDStandardPolicy.</li>

        <li>A user may change some basic behavior of the Marker class by
        providing a special Traits.</li>
    </ul>
 */

#ifndef MARKER_H
#define MARKER_H 1

#include <iostream>
#include <lifev/core/LifeV.hpp>
namespace LifeV
{

//! markerID_Type is the type used to store the geometric entity marker IDs
/*!
 *  An entity marker ID is an integral type that is used to store information about each geometric entity
 *  belonging to a mesh. In particular, it is the number given by the mesh generator that has
 *  generated the mesh to identify different portion of the boundary. It must be convertible to
 *  an ID type.
*/
typedef ID markerID_Type;

//! MarkerIDStandardPolicy - Class that defines the standard policies on Marker Ids
/*!
    This class defines NULLFLAG and how to handle ambiguities among markerIDs
    In particular what to do if a geometric item has to inherit its Marker ID by adjacent items and
    the marker Ids are different. The policy is passed as template argument to the Marker class.
 */
class MarkerIDStandardPolicy
{
public:

    /*! It is the value indicating a null marker ID, i.e am ID not yet
      set to any usable value
    */
    static const markerID_Type S_NULLMARKERID;

    //! Selects the stronger Marker ID
    /*!
     * A dimensional geometric entity G_i may inherit the stronger Marker ID
     * among adiacent geometric entities of greater dimensions.  For
     * example a boundary point Point with an unset Flag may inherit the
     * strongerst Flag of the adjacent boundary faces.
     * It returns a Null Marker ID if any of the entity a or b ha a null Marker ID
    */
    static markerID_Type strongerMarkerID ( markerID_Type const& a, markerID_Type const& b );

    //! Selects the weaker Marker ID between marker IDs
    /*
      ! A lower dimensional geometric entity G_i may inherit the weaker
      Marker ID among the Vertices of the entity.  For example a boundary
      Face with an unset Marker ID may inherit the weakest Marker ID of its
      Vertices. This method may also be used to attribute the Flag to
      Nodes generated on high order elements.
      It returns a null Marker ID if any of the entity a or b has a Null Marker ID.
    */
    static markerID_Type weakerMarkerID ( markerID_Type const& a, markerID_Type const& b );

    //! Equality operator.
    /*
        It is needed in order to select markers with the same Marker ID
    */
    static bool equalMarkerID (const markerID_Type& a, const markerID_Type& b);

};

//! Marker - Base marker class.
/*!
  It stores an object of markerID_Type which may be used for marking a geometric entity.
  The typical use is to specify boundary conditions or material properties associated with
  the entity. The actual boundary conditions will be handled in the DOF
  class. During the creation of a field, the markers are processed to furnish the
  correct boundary conditions.
  The template argument MarkerIDPolicy defaults to MarkerIDStandardPolicy and it defines the way
  ambiguities in the Marker ID definition are treated.
  The template argument FlagPolicy defaults to EntityFlagStandardPolicy and it defines the way
  ambiguities in the flag definition are treated.

  Marker is a concrete base class which also implements basic tool to operate on the markerID.
  All geometric entities stored in a mesh derives from it, thus
  it may be used to extend the capabilities of a geometric entity, for any purpose.

  The default markers used by RegionMesh are contained in MarkerDefinitions.hpp

 */

template <typename MarkerIDPolicy = MarkerIDStandardPolicy>
class Marker
{
public:

    //! @name Public Types
    //@{
    //! The policy used by this marker
    typedef MarkerIDPolicy MarkerIDPolicy_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    Marker();

    //! Constructor given the marker ID
    explicit Marker ( markerID_Type& p );

    //! Copy Constructor
    Marker ( Marker<MarkerIDPolicy> const& markerBase );

    //! Destructor
    virtual ~Marker()
    {
    }

    //@}

    //! @name Methods
    //@{

    //! It enquires if marker ID is different than the null marker ID
    inline bool isMarkerSet() const;

    //! It enquires if marker ID is different than the null marker ID
    inline bool isMarkerUnset() const;

    //! Compares marker IDs
    /*!
        It returns true if the marker ID is equal to the argument
    */
    inline bool hasEqualMarkerID (markerID_Type const& markerID ) const;

    //! Display information about the marker object
    virtual void showMe ( std::ostream& output = std::cout ) const;

    //@}

    //! @name Setters
    //@{

    //! Set marker to the given value
    inline markerID_Type setMarkerID ( markerID_Type const& markerID );

    //! Set marker to the given value only if unset
    markerID_Type updateMarkerID ( markerID_Type const& markerID );

    //! Sets the marker ID to the stronger marker ID of two given markers
    markerID_Type setStrongerMarkerID ( markerID_Type const& markerID1, markerID_Type const& markerID2 );

    //! Sets the marker ID to the weaker marker ID of two given markers
    markerID_Type setWeakerMarkerID ( markerID_Type const& markerID1, markerID_Type const& markerID2 );

    //! Sets to the strongest marker ID
    /*!
        If marker ID is not set, is sets it to that of the argument, otherwise
        it sets it to  the stronger ID between the stored one
        and the one provided by the argument.
     */
    markerID_Type setStrongerMarkerID ( markerID_Type const& markerID );

    //! Sets to the weaker marker ID
    /*!
        If marker ID is not set, is sets it to that of the argument, otherwise
        it sets it to  the weaker ID between the stored one
        and the one provided by the argument.
     */
    markerID_Type setWeakerMarkerID ( markerID_Type const& markerID );

    //! Put marker to NULLFLAG
    inline void unsetMarkerID(); //const;

    //@}

    //! @name Get Methods
    //@{

    //! Extracts the entityFlag associated to the marked entity
    /*!
     * For historical reason this method is called marker, while it should be called
     * markerID() instead. Refactoring however would involve changing too many files and it has been
     * postponed. Just remember that markerID() does not return a Marker!
     */
    inline markerID_Type markerID() const;

    //! Returns the null marker ID

    static markerID_Type const& nullMarkerID();

    //@}

protected:
    markerID_Type M_markerID;
};



//  ***********************************************************************************************************
//                                           IMPLEMENTATION
//  ***********************************************************************************************************

template <typename MarkerIDPolicy>
Marker<MarkerIDPolicy>::Marker() : M_markerID ( MarkerIDPolicy::S_NULLMARKERID )
{
    // nothing to be done here
}

template <typename MarkerIDPolicy>
Marker<MarkerIDPolicy>::Marker ( markerID_Type& markerID ) : M_markerID ( markerID )
{
    // nothing to be done here
}

template <typename MarkerIDPolicy>
Marker<MarkerIDPolicy>::Marker ( Marker<MarkerIDPolicy> const& markerBase ) : M_markerID ( markerBase.markerID() )
{
    // nothing to be done here
}

template <typename MarkerIDPolicy>
markerID_Type Marker<MarkerIDPolicy>::markerID() const
{
    return M_markerID;
}

template <typename MarkerIDPolicy>
markerID_Type const& Marker<MarkerIDPolicy>::nullMarkerID()
{
    return MarkerIDPolicy::S_NULLMARKERID;
}

template <typename MarkerIDPolicy>
markerID_Type Marker<MarkerIDPolicy>::setMarkerID ( markerID_Type const& markerID )
{
    return M_markerID = markerID;
}

template <typename MarkerIDPolicy>
markerID_Type Marker<MarkerIDPolicy>::updateMarkerID ( markerID_Type const& markerID )
{
    if ( M_markerID == nullMarkerID() )
    {
        return setMarkerID ( markerID );
    }
}

template <typename MarkerIDPolicy>
markerID_Type Marker<MarkerIDPolicy>::setStrongerMarkerID ( markerID_Type const& markerID1,
                                                            markerID_Type const& markerID2 )
{
    return setMarkerID ( MarkerIDPolicy::strongerMarkerID ( markerID1, markerID2 ) );
}

template <typename MarkerIDPolicy>
markerID_Type Marker<MarkerIDPolicy>::setWeakerMarkerID ( markerID_Type const& markerID1,
                                                          markerID_Type const& markerID2 )
{
    return setMarkerID ( MarkerIDPolicy::weakerMarkerID ( markerID1, markerID2 ) );
}

template <typename MarkerIDPolicy>
markerID_Type Marker<MarkerIDPolicy>::setStrongerMarkerID ( markerID_Type const& markerID )
{
    if ( isMarkerUnset() )
    {
        return M_markerID = markerID;
    }
    return setMarkerID ( MarkerIDPolicy::strongerMarkerID ( this->markerID(), markerID ) );
}

template <typename MarkerIDPolicy>
markerID_Type Marker<MarkerIDPolicy>::setWeakerMarkerID ( markerID_Type const& markerID )
{
    if ( isMarkerUnset() )
    {
        return M_markerID = markerID;
    }
    return setMarkerID ( MarkerIDPolicy::weakerMarkerID ( this->markerID(), markerID ) );
}

template <typename MarkerIDPolicy>
bool Marker<MarkerIDPolicy>::isMarkerSet() const
{
    return M_markerID != nullMarkerID();
}

template <typename MarkerIDPolicy>
bool Marker<MarkerIDPolicy>::isMarkerUnset() const
{
    return M_markerID == nullMarkerID();
}

template <typename MarkerIDPolicy>
void Marker<MarkerIDPolicy>::unsetMarkerID()
{
    M_markerID = nullMarkerID();
}

template <typename MarkerIDPolicy>
bool Marker<MarkerIDPolicy>::hasEqualMarkerID (markerID_Type const& markerID) const
{
    return MarkerIDPolicy::EqualFlags (markerID, M_markerID);
}

template <typename MarkerIDPolicy>
void Marker<MarkerIDPolicy>::showMe ( std::ostream& output) const
{
    if ( M_markerID == nullMarkerID() )
    {
        output << "UNSET";
    }
    else
    {
        output << M_markerID;
    }
}

} //Namespace LifeV

#endif /* MARKER_H */
