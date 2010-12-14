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
        class (maybe by inheriting from the Marker_Base class, which
        provides the same interfaces of Marker_Base + the user defined
        stuff.</li>
    </ul>

    The flag in the base class+some utilities to select between two flags
    is provided by traits. In particular,

    <ul>
        <li> The MarkerTraits_Base define the basic (compulsory) interface of
        any user defined Marker class.</li>

        <li> The Marker is a class template whose template argument is a
        MarkerTrait, defaulted to MarkerTraits_Base.</li>

        <li>A user may change some basic behavior of the Marker class by
        providing a special Traits.</li>
    </ul>
 */

#ifndef MARKERS_BASE_H
#define MARKERS_BASE_H 1

#include <iostream>
#include <life/lifecore/life.hpp>

namespace LifeV
{

//! MarkerTraits_Base - Class that defines the basic entities and functions.
/*!
    @todo Change the name of the class in MarkerTraitsBase
 */

class MarkerTraits_Base
{
public:

    //! @name Public Types
    //@{

    //! entityFlag_Type is the type used to store the geometric entity flags
    typedef ID entityFlag_Type;
    // Old typedef to delete
    typedef ID EntityFlag;

    //@}

    /*! Nullflag is the value indicating a null flag, i.e a flag not yet
      set to a usable value
    */
    static const entityFlag_Type NULLFLAG;


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

//! Marker_Base - Base marker class.
/*!
  It stores an integral type, aliased to entityFlag_Type, which may be used
  for marking a geometric entity. The typical use is to specify
  boundary conditions or material properties associated with the entity.
  The actual boundary conditions will be handled in the Dof
  class. During the creation of a field, the markers are post-processed
  to furnish the correct boundary conditions.
  The main interfaces are passed trough the MarkerTraits template argument,
  which is defaulted to MarkerTraits_Base

  @todo Change the name of the class in MarkerBase
 */

template <typename MarkerTraits>
class Marker_Base
{
public:

    //! @name Public Types
    //@{
    typedef typename MarkerTraits::entityFlag_Type entityFlag_Type;
    // Old typedef to be removed
    typedef typename MarkerTraits::EntityFlag EntityFlag;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    Marker_Base();

    //! Constructor given the flag
    explicit Marker_Base( entityFlag_Type & p );

    //! Copy Constructor
    Marker_Base( Marker_Base<MarkerTraits> const & markerBase );

    //! Destructor
    virtual ~Marker_Base()
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

    //! Helper function that prints a marker Flag
    std::ostream & __attribute__ ((__deprecated__)) printFlag( entityFlag_Type const f, std::ostream & output ) const;

    //! Helper function that prints "this" marker flag
    std::ostream & __attribute__ ((__deprecated__)) printFlag( std::ostream & output ) const;

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

    //! Put marker to nullflag
    inline void __attribute__ ((__deprecated__)) markerUnset() const;

    //@}

    //! @name Get Methods
    //@{

    //! Extract marker flag
    inline entityFlag_Type marker() const;

    //! Returns the null flag (cannot be modified)
    inline entityFlag_Type const & nullFlag() const;

    //@}

protected:
    entityFlag_Type M_flag;
};

//! MarkerCommon - Contains only typedefs

template
<class MT>
class MarkerCommon
{
public:

    //! @name Public Types
    //@{

    typedef MT markerTraits_Type;
    typedef Marker_Base<MT> pointMarker_Type;
    typedef Marker_Base<MT> edgeMarker_Type;
    typedef Marker_Base<MT> faceMarker_Type;
    typedef Marker_Base<MT> volumeMarker_Type;
    typedef Marker_Base<MT> regionMarker_Type;

    // Old typedefs to delete
    typedef MT MarkerTraits;
    typedef Marker_Base<MT> PointMarker;
    typedef Marker_Base<MT> EdgeMarker;
    typedef Marker_Base<MT> FaceMarker;
    typedef Marker_Base<MT> VolumeMarker;
    typedef Marker_Base<MT> RegionMarker;

    //@}
};


//  ***********************************************************************************************************
//                                           IMPLEMENTATION
//  ***********************************************************************************************************

template <typename MarkerTraits>
Marker_Base<MarkerTraits>::Marker_Base() : M_flag( MarkerTraits::NULLFLAG )
{
    // nothing to be done here
}

template <typename MarkerTraits>
Marker_Base<MarkerTraits>::Marker_Base( entityFlag_Type & flag ) : M_flag( flag )
{
    // nothing to be done here
}

template <typename MarkerTraits>
Marker_Base<MarkerTraits>::Marker_Base( Marker_Base<MarkerTraits> const & markerBase ) : M_flag( markerBase.marker() )
{
    // nothing to be done here
}

template <typename MarkerTraits>
typename MarkerTraits::entityFlag_Type Marker_Base<MarkerTraits>::marker() const
{
    return M_flag;
}

template <typename MarkerTraits>
typename MarkerTraits::entityFlag_Type const & Marker_Base<MarkerTraits>::nullFlag() const
{
    return MarkerTraits::NULLFLAG;
}

template <typename MarkerTraits>
typename MarkerTraits::entityFlag_Type Marker_Base<MarkerTraits>::setMarker( entityFlag_Type const & flag )
{
    return M_flag = flag;
}

template <typename MarkerTraits>
typename MarkerTraits::entityFlag_Type Marker_Base<MarkerTraits>::updateMarker( entityFlag_Type const & flag )
{
    if ( M_flag == nullFlag() )
        return setMarker( flag );
}

template <typename MarkerTraits>
typename MarkerTraits::entityFlag_Type Marker_Base<MarkerTraits>::setStrongerMarker( entityFlag_Type const & flag1, entityFlag_Type const & flag2 )
{
    return setMarker( MarkerTraits::strongerFlag( flag1, flag2 ) );
}

template <typename MarkerTraits>
typename MarkerTraits::entityFlag_Type Marker_Base<MarkerTraits>::setWeakerMarker( entityFlag_Type const & flag1, entityFlag_Type const & flag2 )
{
    return setMarker( MarkerTraits::weakerFlag( flag1, flag2 ) );
}

template <typename MarkerTraits>
typename MarkerTraits::entityFlag_Type Marker_Base<MarkerTraits>::setStrongerMarker( entityFlag_Type const & flag )
{
    if ( isMarkerUnset() )
        return M_flag = flag;
    return setMarker( MarkerTraits::strongerFlag( this->marker(), flag ) );
}

template <typename MarkerTraits>
typename MarkerTraits::entityFlag_Type Marker_Base<MarkerTraits>::setWeakerMarker( entityFlag_Type const & flag )
{
    if ( isMarkerUnset() )
        return M_flag = flag;
    return setMarker( MarkerTraits::weakerFlag( this->marker(), flag ) );
}

template <typename MarkerTraits>
bool Marker_Base<MarkerTraits>::isMarkerSet() const
{
    return M_flag != nullFlag();
}

template <typename MarkerTraits>
bool Marker_Base<MarkerTraits>::isMarkerUnset() const
{
    return M_flag == nullFlag();
}

template <typename MarkerTraits>
void Marker_Base<MarkerTraits>::unsetMarker()
{
    M_flag = nullFlag();
}

template <typename MarkerTraits>
void Marker_Base<MarkerTraits>::markerUnset() const
{
    M_flag = nullFlag();
}


template <typename MarkerTraits>
bool Marker_Base<MarkerTraits>::hasEqualEntityFlag(entityFlag_Type const & flag) const
{
    return MarkerTraits::EqualFlags(flag,M_flag);
}

template <typename MarkerTraits>
std::ostream & Marker_Base<MarkerTraits>::printFlag( entityFlag_Type const f, std::ostream & output ) const
{
    if ( f == nullFlag() )
        output << "UNSET";
    else
        output << f;
    return output;
}

template <typename MarkerTraits>
std::ostream & Marker_Base<MarkerTraits>::printFlag( std::ostream & output ) const
{
    showMe( output );
    return output;
}

template <typename MarkerTraits>
void Marker_Base<MarkerTraits>::showMe( std::ostream & output) const
{
	if ( M_flag == nullFlag() )
		output << "UNSET";
	else
		output << M_flag;
}

} //Namespace LifeV

#endif /* MARKERS_BASE_H */
