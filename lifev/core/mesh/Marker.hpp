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
        EntityFLAG). The main purpose of the Entity flag is to associate
        boundary conditions or material properties to the Geometry
        entity. </li>

        <li> There is a mechanism, via the MarkerCommon class template, to
        extend the mesh entities definitions by adding user defined data
        and methods. It this case the user should create its own marker
        class (maybe by inheriting from the Marker class, which
        provides the same interfaces of Marker + the user defined
        stuff.</li>
    </ul>

    The flag in the base class+some utilities to select between two flags
    is provided by traits. In particular,

    <ul>
        <li> The EntityFlagStandardPolicy define the basic (compulsory) interface of
        any user defined Marker class.</li>

        <li> The Marker is a class template whose template argument is a
        MarkerTrait, defaulted to EntityFlagStandardPolicy.</li>

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

//! markerID_Type is the type used to store the geometric entity flags
/*!
 *  An entity flag is an integral type that is used to store information about each geometric entity
 *  belonging to a mesh. In particular, it is the number given by the mesh generator that has
 *  generated the mesh to identify different portion of the boundary. It mast be convertible with
 *  an ID type.
*/
typedef ID markerID_Type;

//! EntityFlagStandardPolicy - Class that defines the standard policies on EntityFlags
/*!
    This class defines NULLFLAG and how to handle ambiguities among markerIDs
    In particular what to do if a geometric item has to inherit its flag by adjacent items and
    the flags are different. The policy is passed as template argument to the Marker class.
 */
class EntityFlagStandardPolicy
{
public:

    /*! Nullflag is the value indicating a null flag, i.e a flag not yet
      set to any usable value
    */
    static const markerID_Type S_NULLFLAG;

    //! Selects the stronger between two flags
    /*! A dimensional geometric entity G_i may inherit the stronger flag
      among adjacent geometric entities of greater dimensions.  For
      example a boundary point Point with an unset Flag may inherit the
      strongest Flag of the adjacent boundary faces.
      It returns NULLFLAG if any of the entity a or b is a NULLFLAG.
    */
    static markerID_Type strongerFlag( markerID_Type const & a, markerID_Type const & b );

    //! Selects the weaker between two flags
    /*
      ! A lower dimensional geometric entity G_i may inherit the weaker
      flag among the Vertices of the entity.  For example a boundary
      Face with an unset Flag may inherit the weakest Flag of its
      Vertices. This method may also be used to attribute the Flag to
      Nodes generated on high order elements.
      It returns NULLFLAG if any of the entity a or b is a NULLFLAG.
    */
    static markerID_Type weakerFlag( markerID_Type const & a, markerID_Type const & b );

    //! Equality operator.
    /*
        It is needed in order to select markers with the same entity flag
    */
    static bool EqualFlags(const markerID_Type& a, const markerID_Type& b);

};

//! Marker - Base marker class.
/*!
  It stores an object of markerID_Type which may be used for marking a geometric entity.
  The typical use is to specify boundary conditions or material properties associated with
  the entity. The actual boundary conditions will be handled in the DOF
  class. During the creation of a field, the markers are processed to furnish the
  correct boundary conditions.
  The template argument FlagPolicy defaults to EntityFlagStandardPolicy and it defines the way
  ambiguities in the flag definition are treated.

  Marker is a concrete base class which also implements basic tool to operate on the markerID.
  All geometric entities stored in a mesh derives from it, thus
  it may be used to extend the capabilities of a geometric entity, for any purpose.

  The default markers used by RegionMesh are contained in MarkerDefinitions.hpp

 */

template <typename FlagPolicy=EntityFlagStandardPolicy>
class Marker
{
public:

    //! @name Public Types
    //@{
    //! The policy used by this marker
    typedef FlagPolicy flagPolicy_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    Marker();

    //! Constructor given the flag
    explicit Marker( markerID_Type & p );

    //! Copy Constructor
    Marker( Marker<FlagPolicy> const & markerBase );

    //! Destructor
    virtual ~Marker()
    {
    }

    //@}

    //! @name Methods
    //@{

    //! It enquires if marker flag is different than the nullflag
    inline bool isMarkerSet() const;

    //! It enquires if marker flag is different than the nullflag
    inline bool isMarkerUnset() const;

    //! Compares flags
    /*!
        It returns true if the marker flag is equal to the argument
    */
    inline bool hasEqualEntityFlag(markerID_Type const & flag ) const;

    //! Display information about the marker object
    virtual void showMe( std::ostream & output = std::cout ) const;

    //@}

    //! @name Setters
    //@{

    //! Set marker to the given value
    inline markerID_Type setMarker( markerID_Type const & flag );

    //! Set marker to the given value only if unset
    markerID_Type updateMarker( markerID_Type const & flag );

    //! Sets the flag to the stronger flag of two given markers
    markerID_Type setStrongerMarker( markerID_Type const & flag1, markerID_Type const & flag2 );

    //! Sets the flag to the weaker flag of two given markers
    markerID_Type setWeakerMarker( markerID_Type const & flag1, markerID_Type const & flag2 );

    //! Sets to the strongest flag
    /*!
        If marker flag is not set, is sets it to that of the argument, otherwise
        it sets it to  the stronger flag between the stored one
        and the one provided by the argument.
     */
    markerID_Type setStrongerMarker( markerID_Type const & flag );

    //! Sets to the strongest flag
    /*!
        If marker flag is not set, is sets it to that of the argument, otherwise
        it sets it to  the weaker flag between the stored one
        and the one provided by the argument.
     */
    markerID_Type setWeakerMarker( markerID_Type const & flag );

    //! Put marker to NULLFLAG
    inline void unsetMarker(); //const;

    //@}

    //! @name Get Methods
    //@{

    //! Extracts the entityFlag associated to the marked entity
    /*!
     * For historical reason this method is called marker, while it should be called
     * markerID() instead. Refactoring however would involve changing too many files and it has been
     * postponed. Just remember that marker() does not return a Marker!
     */
    inline markerID_Type marker() const;

    //! Returns the null flag

    static markerID_Type const & nullFlag();

    //@}

protected:
    markerID_Type M_flag;
};



//  ***********************************************************************************************************
//                                           IMPLEMENTATION
//  ***********************************************************************************************************

template <typename FlagPolicy>
Marker<FlagPolicy>::Marker() : M_flag( FlagPolicy::S_NULLFLAG )
{
    // nothing to be done here
}

template <typename FlagPolicy>
Marker<FlagPolicy>::Marker( markerID_Type & flag ) : M_flag( flag )
{
    // nothing to be done here
}

template <typename FlagPolicy>
Marker<FlagPolicy>::Marker( Marker<FlagPolicy> const & markerBase ) : M_flag( markerBase.marker() )
{
    // nothing to be done here
}

template <typename FlagPolicy>
markerID_Type Marker<FlagPolicy>::marker() const
{
    return M_flag;
}

template <typename FlagPolicy>
markerID_Type const & Marker<FlagPolicy>::nullFlag()
{
    return FlagPolicy::S_NULLFLAG;
}

template <typename FlagPolicy>
markerID_Type Marker<FlagPolicy>::setMarker( markerID_Type const & flag )
{
    return M_flag = flag;
}

template <typename FlagPolicy>
markerID_Type Marker<FlagPolicy>::updateMarker( markerID_Type const & flag )
{
    if ( M_flag == nullFlag() )
        return setMarker( flag );
}

template <typename FlagPolicy>
markerID_Type Marker<FlagPolicy>::setStrongerMarker( markerID_Type const & flag1, markerID_Type const & flag2 )
{
    return setMarker( FlagPolicy::strongerFlag( flag1, flag2 ) );
}

template <typename FlagPolicy>
markerID_Type Marker<FlagPolicy>::setWeakerMarker( markerID_Type const & flag1, markerID_Type const & flag2 )
{
    return setMarker( FlagPolicy::weakerFlag( flag1, flag2 ) );
}

template <typename FlagPolicy>
markerID_Type Marker<FlagPolicy>::setStrongerMarker( markerID_Type const & flag )
{
    if ( isMarkerUnset() )
        return M_flag = flag;
    return setMarker( FlagPolicy::strongerFlag( this->marker(), flag ) );
}

template <typename FlagPolicy>
markerID_Type Marker<FlagPolicy>::setWeakerMarker( markerID_Type const & flag )
{
    if ( isMarkerUnset() )
        return M_flag = flag;
    return setMarker( FlagPolicy::weakerFlag( this->marker(), flag ) );
}

template <typename FlagPolicy>
bool Marker<FlagPolicy>::isMarkerSet() const
{
    return M_flag != nullFlag();
}

template <typename FlagPolicy>
bool Marker<FlagPolicy>::isMarkerUnset() const
{
    return M_flag == nullFlag();
}

template <typename FlagPolicy>
void Marker<FlagPolicy>::unsetMarker()
{
    M_flag = nullFlag();
}

template <typename FlagPolicy>
bool Marker<FlagPolicy>::hasEqualEntityFlag(markerID_Type const & flag) const
{
    return FlagPolicy::EqualFlags(flag,M_flag);
}

template <typename FlagPolicy>
void Marker<FlagPolicy>::showMe( std::ostream & output) const
{
	if ( M_flag == nullFlag() )
		output << "UNSET";
	else
		output << M_flag;
}

} //Namespace LifeV

#endif /* MARKER_H */
