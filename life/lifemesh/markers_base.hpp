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
#ifndef HH_MARKERSBASE_HH_
#define HH_MARKERSBASE_HH_
#include <iostream>

namespace LifeV
{
/*            ************************ MARKER BASE **************************/
/*! \file marker_base.h Basic definition of markers

Here we define the basic markers. Markers have two purposes:

<ul>

  <li> To add an indicator (EntityFlag) to all geometry entities. In the
  base MarkerTrait this indicator is a long integer (aliased to
  EntityFLAG). The main purpose of the Entity flag is to associate
  boundary conditions or material properties to the Geometry
  entity. </li>


  <li> There is a mechanism, via the MarkerCommon class template, to
  extend the mesh enitities definitions by adding user defined data
  and methods. It this case the user should create its own marker
  class (maybe by inheriting from the Marker_Base class, which
  provides the same interfaces of Marker_Base + the user defined
  stuff.  </li>

</ul>

The flag in the base class+some utilities to select between two flags
is provided by traits. In particular, <ul>

<li> The MarkerTraits_Base define the basic (compulsory) interface of
any user defined Marker class.

  <li> The Marker is a class template whose template argument is a
  MarkerTrait, defaulted to MarkerTraits_Base.

<li>A user may change some basic behaviour of the Marker class by providing a special Traits.
</ul>

*/


//! This defines the basic entities and functions. The user may change them
//! if needed

class MarkerTraits_Base
{
public:
  //! EntityFlag is the type used to store the geometric entity flags
  /*!  It is a signed long int so thet we have a lot of room */

  typedef signed long int EntityFlag;

  /*! Nullflag is the value indicating a null flag, i.e a flag not yet
    set to a usable value*/

  static const EntityFlag NULLFLAG;

  //! Selects the stronger between two flags
  /*! A dimensional geometric entity G_i may inherit the stronger flag
    among adiacent geometric entities of greater dimensions.  For
    example a boundary point Point with an unset Flag may inherit the
    strongerst Flag of the adjacent boundary faces.
    It returns NULLFLAG if any of the entity a or b is a NULLFLAG.
  */

  static EntityFlag strongerFlag(EntityFlag const & a,EntityFlag const & b);

  //! Selects the weaker between two flags
  /*
    ! A lower dimensional geometric entity G_i may inherit the weaker
    flag among the Vertices of the entity.  For example a boundary
    Face with an unset Flag may inherit the weakest Flag of its
    Vertices. This method may also be used to attribute the Flag to
    Nodes generated on high order elements.
    It returns NULLFLAG if any of the entity a or b is a NULLFLAG.
  */
  static EntityFlag weakerFlag(EntityFlag const & a,EntityFlag const & b);
};

// const MarkerTraits_Base::EntityFlag MarkerTraits_Base::NULLFLAG=LONG_MIN;




//! Base marker class.
/*!
  It stores an integral type, aliased to EntityFlag, which may be used
  for marking a geometrical entity. The typical use is to specify
  boundary conditions or material properties associated with the entity.
  The actual boundary conditions will be handled in the Dof
  class. During the creation of a field, the markers are post-processed
  to furnish the correct boundary conditions.
  The main interfaces are passed trough the MarkerTraits template argiment, which is defaulted to MarkerTraits_Base*/

template<typename MarkerTraits>
class Marker_Base {
public:

  typedef typename MarkerTraits::EntityFlag EntityFlag;

  //! Null flag as default
  explicit Marker_Base();

  //! Give the flag
  explicit Marker_Base(EntityFlag & p);

  Marker_Base(Marker_Base<MarkerTraits> const & m);

  //! Extract marker flag
  inline EntityFlag marker() const;

  //! Returns the null flag (cannot be modified)
  inline EntityFlag const & nullFlag() const;

  //! Set marker to given value
  inline EntityFlag setMarker(EntityFlag const & c);

  //! Set marker to given value only if unset
  EntityFlag updateMarker(EntityFlag const & c);

  //! Sets the flag to the stronger flag of two given markers

  EntityFlag setStrongerMarker(EntityFlag const & p1, EntityFlag const & p2);

  //! Sets the flag to the weaker flag of two given markers

  EntityFlag setWeakerMarker(EntityFlag const & p1, EntityFlag const & p2);

  //! If marker flag is unset, is stes it to that of the argument, otherwise
  //! is sets it to  the stronger flag between the stored one
  //! and the one provided by the argument.
  EntityFlag setStrongerMarker(EntityFlag const & p);

  //! If marker flag is unset , it sets it to that of the argument, otherwise
  //! is sets it to  the weaker flag between the stored one
  //! and the one provided by the argument.

  EntityFlag setWeakerMarker(EntityFlag const & p);

  //! It enquires if marker flag is different than the nullflag
  inline bool isMarkerSet() const;

  //! It enquires if marker flag is different than the nullflag
  inline bool isMarkerUnset() const;

  //! Put marker to nullflag
  inline bool unsetMarker() const;

  //! Put marker to nullflag
  inline bool markerUnset() const;

  //! Helper function that prints a marker Flag
  std::ostream & printFlag(EntityFlag const f, std::ostream & out) const;

  //! Helper function that prints "this" marker flag
  std::ostream &  printFlag(std::ostream & out) const;

protected:
  EntityFlag flag;
};

//! Defaults markers
// typedef Marker_Base<MarkerTraits_Base> DefPointMarker;
// typedef Marker_Base<MarkerTraits_Base> DefEdgeMarker;
// typedef Marker_Base<MarkerTraits_Base> DefFaceMarker;
// typedef Marker_Base<MarkerTraits_Base> DefVolumeMarker;
// typedef Marker_Base<MarkerTraits_Base> DefRegionMarker;

//! Class containing the markers
// template
// <
// class PM=DefPointMarker,
//   class EM=DefEdgeMarker,
//   class FM=DefFaceMarker,
//   class VM=DefVolumeMarker,
//   class RM=DefRegionMarker
// >
// class MarkerCommon
// {
// public:
//   typedef PM PointMarker;
//   typedef EM EdgeMarker;
//   typedef FM FaceMarker;
//   typedef VM VolumeMarker;
//   typedef RM RegionMarker;
// };

template
<class MT>
class MarkerCommon
{
public:
  typedef MT MarkerTraits;
  typedef Marker_Base<MT> PointMarker;
  typedef Marker_Base<MT> EdgeMarker;
  typedef Marker_Base<MT> FaceMarker;
  typedef Marker_Base<MT> VolumeMarker;
  typedef Marker_Base<MT> RegionMarker;
};


//  ***********************************************************************************************************
//                                           IMPLEMENTATION
//  ***********************************************************************************************************
//Marker_Base<MarkerTraits>

template<typename MarkerTraits>
Marker_Base<MarkerTraits>::Marker_Base():flag(MarkerTraits::NULLFLAG)
{
    // nothing to be done here
}

template<typename MarkerTraits>
Marker_Base<MarkerTraits>::Marker_Base(EntityFlag & p):flag(p)
{
    // nothing to be done here
}

template<typename MarkerTraits>
Marker_Base<MarkerTraits>::Marker_Base(Marker_Base<MarkerTraits> const & m):flag(m.marker())
{
    // nothing to be done here
}

template<typename MarkerTraits>
typename MarkerTraits::EntityFlag Marker_Base<MarkerTraits>::marker() const
{
    return flag;
}

template<typename MarkerTraits>
typename MarkerTraits::EntityFlag const & Marker_Base<MarkerTraits>::nullFlag() const
{
    return MarkerTraits::NULLFLAG;
}

template<typename MarkerTraits>
typename MarkerTraits::EntityFlag Marker_Base<MarkerTraits>::setMarker(EntityFlag const & c){return flag=c;}

template<typename MarkerTraits>
typename MarkerTraits::EntityFlag Marker_Base<MarkerTraits>::updateMarker(EntityFlag const & c){if (flag==nullFlag())
  return setMarker(c); }

template<typename MarkerTraits>
typename MarkerTraits::EntityFlag Marker_Base<MarkerTraits>::setStrongerMarker(EntityFlag const & p1, EntityFlag const & p2)
{
  return setMarker(MarkerTraits::strongerFlag(p1,p2));
}

template<typename MarkerTraits>
typename MarkerTraits::EntityFlag Marker_Base<MarkerTraits>::setWeakerMarker(EntityFlag const & p1, EntityFlag const & p2)
{
  return setMarker(MarkerTraits::weakerFlag(p1,p2));
}

template<typename MarkerTraits>
typename MarkerTraits::EntityFlag Marker_Base<MarkerTraits>::setStrongerMarker(EntityFlag const & p)
{
  if (markerUnset()) return flag=p;
  return setMarker(MarkerTraits::strongerFlag(this->marker(),p));
}

template<typename MarkerTraits>
typename MarkerTraits::EntityFlag Marker_Base<MarkerTraits>::setWeakerMarker(EntityFlag const & p)
{
  if (markerUnset()) return flag=p;
  return setMarker(MarkerTraits::weakerFlag(this->marker(),p));
}

template<typename MarkerTraits>
bool Marker_Base<MarkerTraits>::isMarkerSet() const {return flag != nullFlag();}

template<typename MarkerTraits>
bool Marker_Base<MarkerTraits>::isMarkerUnset() const {return flag == nullFlag();}

template<typename MarkerTraits>
bool Marker_Base<MarkerTraits>::unsetMarker() const {return isMarkerUnset();}


template<typename MarkerTraits>
bool Marker_Base<MarkerTraits>::markerUnset() const {return isMarkerUnset();}

template<typename MarkerTraits>
std::ostream & Marker_Base<MarkerTraits>::printFlag(EntityFlag const f, std::ostream & out) const{
  if (f ==nullFlag()) out <<"UNSET"; else out << f;
  return out;
}

template<typename MarkerTraits>
std::ostream & Marker_Base<MarkerTraits>::printFlag(std::ostream & out) const{
  return printFlag(flag,out);
}
}

#endif

