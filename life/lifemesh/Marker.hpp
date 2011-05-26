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
        <li> To add an indicator (entityFlag_Type) to all geometry entities. In the
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
#include <life/lifecore/LifeV.hpp>

namespace LifeV
{

//! entityFlag_Type is the type used to store the geometric entity flags
/*!
 *  An entity flag is an integral type that is used to store information about each geometric entity
 *  belonging to a mesh. In particular, it is the number given by the mesh generator that has
 *  generated the mesh to identify different portion of the boundary.
*/
typedef ID entityFlag_Type;


//! EntityFlagStandardPolicy - Class that defines the standard policies on EntityFlags
/*!
    This class defines the NULLFLAG and how to handle ambiguities among EntityFlags
    In particular what to do if a geometric item has to inherit its flag by adjacent items and
    the flags are different
 */
class EntityFlagStandardPolicy
{
public:

    //! @name Public Types
    //@{

    //! entityFlag_Type is the type used to store the geometric entity flags
    typedef ID entityFlag_Type;

    //@}

    /*! Nullflag is the value indicating a null flag, i.e a flag not yet
      set to a usable value
    */
    static const entityFlag_Type S_NULLFLAG;

    //! Selects the stronger between two flags
    /*! A dimensional geometric entity G_i may inherit the stronger flag
      among adiacent geometric entities of greater dimensions.  For
      example a boundary point Point with an unset Flag may inherit the
      strongerst Flag of the adjacent boundary faces.
      It returns NULLFLAG if any of the entity a or b is a NULLFLAG.
    */
    static entityFlag_Type strongerFlag( entityFlag_Type const & a, entityFlag_Type const & b );

    //! Selects the weaker between two flags
    /*
      ! A lower dimensional geometric entity G_i may inherit the weaker
      flag among the Vertices of the entity.  For example a boundary
      Face with an unset Flag may inherit the weakest Flag of its
      Vertices. This method may also be used to attribute the Flag to
      Nodes generated on high order elements.
      It returns NULLFLAG if any of the entity a or b is a NULLFLAG.
    */
    static entityFlag_Type weakerFlag( entityFlag_Type const & a, entityFlag_Type const & b );

    //! Equality operator.
    /*
        It is needed in order to select markers with the same entity flag
    */
    static bool EqualFlags(const entityFlag_Type& a, const entityFlag_Type& b);

};

//! Marker - Base marker class.
/*!
  It stores an integral type, of entityFlag_Type (we indicate it by entityflag),
  which may be used for marking a geometric entity. The typical use is to specify
  boundary conditions or material properties associated with the entity.
  The actual boundary conditions will be handled in the Dof
  class. During the creation of a field, the markers are post-processed
  to furnish the correct boundary conditions.
  The main interfaces are passed trough the EntityFlagStandardPolicy template argument,
  which is defaulted to EntityFlagStandardPolicy

  It is a concrete base class which also implements basic tool to operate on the entityflag
  using the policy passed as template argument. Since all geometric entities stored in a mesh
  derives from it, it may be used to extend the capabilities of a goemetric entity, for any purpose.

  The default markers used by RegionMesh are contained in MarkerDefinitions.hpp

 */

template <typename MarkerTraits>
class Marker
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    Marker();

    //! Constructor given the flag
    explicit Marker( entityFlag_Type & p );

    //! Copy Constructor
    Marker( Marker<MarkerTraits> const & markerBase );

    //! Destructor
    virtual ~Marker()
    {
        // Nothing to be done
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
    inline bool hasEqualEntityFlag(entityFlag_Type const & flag ) const;

    //! Display method
    virtual void showMe( std::ostream & output = std::cout ) const;

    //@}

    //! @name Set Methods
    //@{

    //! Set marker to given value
    inline entityFlag_Type setMarker( entityFlag_Type const & flag );

    //! Set marker to given value only if unset
    entityFlag_Type updateMarker( entityFlag_Type const & flag );

    //! Sets the flag to the stronger flag of two given markers
    entityFlag_Type setStrongerMarker( entityFlag_Type const & flag1, entityFlag_Type const & flag2 );

    //! Sets the flag to the weaker flag of two given markers
    entityFlag_Type setWeakerMarker( entityFlag_Type const & flag1, entityFlag_Type const & flag2 );

    //! Sets to the strongest flag
    /*!
        If marker flag is not set, is sets it to that of the argument, otherwise
        it sets it to  the stronger flag between the stored one
        and the one provided by the argument.
     */
    entityFlag_Type setStrongerMarker( entityFlag_Type const & flag );

    //! Sets to the strongest flag
    /*!
        If marker flag is not set, is sets it to that of the argument, otherwise
        it sets it to  the weaker flag between the stored one
        and the one provided by the argument.
     */
    entityFlag_Type setWeakerMarker( entityFlag_Type const & flag );

    //! Put marker to nullflag
    inline void unsetMarker(); //const;

    //@}

    //! @name Get Methods
    //@{

    //! Extracts the enitytFlag associated to the marked entity
    /*!
     * For hystorical reason this method is called marker, while it should be called
     * entityFlag() instead. Refactoring however would involve changing too many files and it has been
     * posponed. Just remember that marker() does not return a Marker!
     */
    inline entityFlag_Type marker() const;

    //! Returns the null flag

    inline entityFlag_Type const & nullFlag() const;

    //@}

protected:
    entityFlag_Type M_flag;
};



//  ***********************************************************************************************************
//                                           IMPLEMENTATION
//  ***********************************************************************************************************

template <typename MarkerTraits>
Marker<MarkerTraits>::Marker() : M_flag( MarkerTraits::S_NULLFLAG )
{
    // nothing to be done here
}

template <typename MarkerTraits>
Marker<MarkerTraits>::Marker( entityFlag_Type & flag ) : M_flag( flag )
{
    // nothing to be done here
}

template <typename MarkerTraits>
Marker<MarkerTraits>::Marker( Marker<MarkerTraits> const & markerBase ) : M_flag( markerBase.marker() )
{
    // nothing to be done here
}

template <typename MarkerTraits>
entityFlag_Type Marker<MarkerTraits>::marker() const
{
    return M_flag;
}

template <typename MarkerTraits>
entityFlag_Type const & Marker<MarkerTraits>::nullFlag() const
{
    return MarkerTraits::S_NULLFLAG;
}

template <typename MarkerTraits>
entityFlag_Type Marker<MarkerTraits>::setMarker( entityFlag_Type const & flag )
{
    return M_flag = flag;
}

template <typename MarkerTraits>
entityFlag_Type Marker<MarkerTraits>::updateMarker( entityFlag_Type const & flag )
{
    if ( M_flag == nullFlag() )
        return setMarker( flag );
}

template <typename MarkerTraits>
entityFlag_Type Marker<MarkerTraits>::setStrongerMarker( entityFlag_Type const & flag1, entityFlag_Type const & flag2 )
{
    return setMarker( MarkerTraits::strongerFlag( flag1, flag2 ) );
}

template <typename MarkerTraits>
entityFlag_Type Marker<MarkerTraits>::setWeakerMarker( entityFlag_Type const & flag1, entityFlag_Type const & flag2 )
{
    return setMarker( MarkerTraits::weakerFlag( flag1, flag2 ) );
}

template <typename MarkerTraits>
entityFlag_Type Marker<MarkerTraits>::setStrongerMarker( entityFlag_Type const & flag )
{
    if ( isMarkerUnset() )
        return M_flag = flag;
    return setMarker( MarkerTraits::strongerFlag( this->marker(), flag ) );
}

template <typename MarkerTraits>
entityFlag_Type Marker<MarkerTraits>::setWeakerMarker( entityFlag_Type const & flag )
{
    if ( isMarkerUnset() )
        return M_flag = flag;
    return setMarker( MarkerTraits::weakerFlag( this->marker(), flag ) );
}

template <typename MarkerTraits>
bool Marker<MarkerTraits>::isMarkerSet() const
{
    return M_flag != nullFlag();
}

template <typename MarkerTraits>
bool Marker<MarkerTraits>::isMarkerUnset() const
{
    return M_flag == nullFlag();
}

template <typename MarkerTraits>
void Marker<MarkerTraits>::unsetMarker()
{
    M_flag = nullFlag();
}

template <typename MarkerTraits>
bool Marker<MarkerTraits>::hasEqualEntityFlag(entityFlag_Type const & flag) const
{
    return MarkerTraits::EqualFlags(flag,M_flag);
}

template <typename MarkerTraits>
void Marker<MarkerTraits>::showMe( std::ostream & output) const
{
	if ( M_flag == nullFlag() )
		output << "UNSET";
	else
		output << M_flag;
}

} //Namespace LifeV

#endif /* MARKER_H */
