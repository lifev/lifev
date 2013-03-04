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
    @brief It contains the standard selector for internal entities

    @author
    @contributor Nur Aiman Fadel <nur.fadel@mail.polimi.it>
    @maintainer Nur Aiman Fadel <nur.fadel@mail.polimi.it>

    @date

    A more detailed description of the file (if necessary)
 */

#ifndef _SELECTMARKER_HH_
#define _SELECTMARKER_HH_ 1

#include <functional>
#include <lifev/core/mesh/Marker.hpp>
#include <lifev/core/mesh/MeshEntityContainer.hpp>

namespace LifeV
{

/** InternalEntitySelector.
 *
 *  Functor class that tells whether a marker ID corresponds to an internal face
 *  @author Luca Formaggia
 *  @see SetFlagAccordingToWatermarks
 *
 *   This class takes in input an EntityFlag and it can be used <br>
 *   in order to understand if the input is or not an internal face.
 *
 *   @deprecated OBSOLETE, now this is handled through SetFlagAccordingToWatermarks
 */
class InternalEntitySelector
{
public:
    //! The default watermark used when standard constructor is adopted.
    static const markerID_Type defMarkFlag;

    /** @name Constructors & Destructor */
    //@{

    //! Empty Constructor
    InternalEntitySelector();

    //! Short description of the constructor
    /*!
        It is a constructor which requires the EntityFlag
        @param w, the costant EntityFlag which is required in order
        to create an InternalEntitySelector object.
     */
    explicit InternalEntitySelector ( const markerID_Type& w );
    //@}

    //! @name Operators
    //@{

    //! The round brackets operator
    /*!
        Operator returning true if the flag corresponds to internal entity.
        If the EntityFlag is greater that the watermark, the associated geometry entity is internal.
        @param test, it is the reference to geometric entity.
        @return true, if the flag corresponds to an internal entity.
     */
    bool operator() ( markerID_Type const& test ) const;
    //@}

private:

    //! The current watermark
    markerID_Type M_watermarkFlag;

}; // Class InternalEntitySelector

namespace Utility
{

//! A helper functor
/*!
 * This is a helper functor not of general use.
 */
class MarkerMapTraits
{
public:
    //! The type of a range of marker ID values
    typedef std::pair<markerID_Type, markerID_Type> rangeID_Type;

    //! Comparison operator for ranges
    bool operator() ( rangeID_Type const& a, rangeID_Type const& b ) const
    {
        return a.first < b.first;
    }
};

}// end namespace Utility

/** @defgroup MeshEntityFlagsChangers
 * Useful functors to change marker IDs according to certain conditions
 *
 * @{
 */
//!Set the flags of a mesh entity according to predefined ranges
/*!
*  This utility allows to associate to ranges [markerID_Type a, markerID_Type b] a marker ID.
*  Then the operator ()(MeshEntitity const & e) will change the flag according to the following rule
*  -# if e.marker falls in a stored range (a<=e.markerID()<=b) the corresponding flag is used
*     with the given policy (which defaults to turnOn, so the flag is activated)
*  -# otherwise the flag is left unchanged
*
*  An example of its usage
*  @verbatim
*  setFlagAccordingToMarkerRanges changer(Flag::turnOn); //I may use the default constructor
*  changer.insert(std::make_pair(3,5),entityFlags::INTERNAL_INTERFACE);
*  changer.insert(std::make_pair(8,10),entityFlags::VERTEX);
*  mesh.pointList.changeAccordingToFunctor(changer);
*  @endverbatim
*
*  Now all points with marker in the range [3,5] are set as INTERNAL_INTERFACE and those in [8,10]
*  to VERTEX.
*
*  @note It is up to the user to give consistent, i.e. non overlapping, ranges.
*
*/

class SetFlagAccordingToMarkerRanges
{
public:
    //! The type of a flag policy
    typedef flag_Type (*flagPolicy_ptr) (flag_Type const&, flag_Type const&);
    //! The type of a range
    typedef Utility::MarkerMapTraits::rangeID_Type rangeID_Type;

    //! Constructor optionally takes a policy
    explicit SetFlagAccordingToMarkerRanges ( const flagPolicy_ptr& flagPolicy = &Flag::turnOn ) :
        M_flagPolicy (flagPolicy) {};

    //! Inserts a range with associated flag
    /*!
     * The range may be given by calling std::make_pair(markerID_Type a, markerID_Type b)
     * It is the user responsibility making sure that b>=a and that there are not overlapping
     * ranges. Otherwise the result is unpredictable
     *
     * @param key a pair of markerIDs defining a range
     * @param flag the associated flag
     */
    void insert (rangeID_Type const& key, flag_Type flag);

    //! Finds the flag associated to a range.
    /*!
     * It is mainly meant as helper function for operator(), but it is exposed
     * because it might be useful for debugging
     * @param m The markerID to be searched
     * @return a pair containing the flag and a bool. If the bool is false it means that there is no
     *         corresponding range. In that case the flag defaults to DEFAULT
     */
    std::pair<flag_Type, bool> findFlag ( markerID_Type const& m) const;

    //! The operator doing the job
    /*!
     * It is the operator that does the check on a given mesh entity
     * @param e The entity to change (possibly) using the given policy
     */
    template<typename MeshEntity>
    void operator() (MeshEntity& e) const
    {
        std::pair<flag_Type, bool> tmp = this->findFlag ( e.markerID() );
        if ( tmp.second )
        {
            e.replaceFlag ( M_flagPolicy ( e.flag(), tmp.first ) );
        }
    }

private:
    typedef std::map<rangeID_Type, flag_Type, Utility::MarkerMapTraits> map_Type;
    typedef map_Type::const_iterator const_iterator_Type;
    map_Type M_map;
    flagPolicy_ptr M_flagPolicy;
};

//!Sets the flag of a mesh entity according to the value of its Marker id.
/*!
*  We compare the marker ID with a watermark according to a policy which is a comparison operator
*
*   @code
    bool operator()(markerID_Type const & mId, markerID_Type const & watermark)
*   @endcode
*
*  which by default is greater<>  (i.e. returns true if  mId > watermark).
*  If the comparison is true the given flag is assigned to the entity using a second policy
*  implemented as a function:
*
*  @code
   flag_Type flagPolicy  ( flag_Type const & inputFlag, flag_Type const & refFlag )
*  @endcode
*
*  passed in the constructor and defaulted to Flag::turnOn
*
*  Example:
*  I want to set the flag INTERNAL_INTERFACE to all faces with marker ID > 1000
*
* @code
   SetFlagAccordingToWatermark
          changer(EntityFlags::INTERNAL_INTERFACE,1000,Flag::turnOn)
   // the last argument (turnOn) is not needed
   mesh.faceList.changeAccordingToFunctor(changer);
  @endcode
*  I want to set the flag INTERNAL_INTERFACE to all faces with marker ID =12000
*  @code
   SetFlagAccordingToWatermark<std::equal_to<markerID_Type> >
                           changer(EntityFlags::INTERNAL_INTERFACE,12000,Flag::turnOn)
   (the last argument is not needed)
   mesh.faceList.changeAccordingToFunctor(changer);
*  @endcode
*
*/
template<typename Policy = std::greater<markerID_Type> >
class SetFlagAccordingToWatermark
{
public:
    typedef flag_Type (*flagPolicy_ptr) ( flag_Type const&, flag_Type const& );

    //! The constructor
    /*!
     * @param flagToSet the flag we want to set if the condition is satisfied
     * @param watermark The watermark to be used in the comparison
     * @param flagPolicy The policy we want to use. Default is turnOn, which
     * turns on the bit-flags defined in flagToSet. Other possibilities are turnOff or
     * change (the names are self explanatory)
     */
    SetFlagAccordingToWatermark ( const flag_Type& flagToSet,
                                  const markerID_Type& watermark,
                                  const flagPolicy_ptr& flagPolicy = &Flag::turnOn ) :
        M_flagToSet ( flagToSet ),
        M_watermark ( watermark ),
        M_policy ( Policy() ),
        M_flagPolicy ( flagPolicy ) {}

    //! The operator doing the job
    /*!
     * It is the operator that does the check on a given mesh entity
     * @param e The entity to change (possibly) using the given policy
     */
    template<typename MeshEntity>
    void operator() ( MeshEntity& e ) const
    {
        if ( M_policy (e.markerID(), M_watermark) )
        {
            e.replaceFlag ( M_flagPolicy ( e.flag(), M_flagToSet ) );
        }
    }

private:
    const flag_Type M_flagToSet;
    const markerID_Type M_watermark;
    const Policy M_policy;
    const flagPolicy_ptr M_flagPolicy;
};

//! Sets the flag of a mesh entity according to the value of a vector of Marker IDs.
/*!
 *  We compare the marker IDs with all values contained in the vector using the
 *  equal_to<> operator.
 *  If the comparison is true the given flag is assigned to the entity using the
 *  policy implemented as a function
 *  @code
    flag_Type flagPolicy  ( flag_Type const & inputFlag, flag_Type const & refFlag )
 *  @endcode
 *  passed in the constructor and defaulted to Flag::turnOn
 *
 *  Example:
 *  I want to set the flag INTERNAL_INTERFACE to all faces with marker ID = 1000, 2000 and 3000
 * @code
   vector<markerID_Type> fl; fl.push_back(1000); fl.push_back(2000); fl.push_back(3000)
   SetFlagAccordingToWatermarks<face_Type> changer(INTERNAL_INTERFACE,fl,Flag::turnOn)
   //the last argument is not needed
   mesh.faceList.changeAccordingToFunctor(changer);
   @code
 *
 */
class SetFlagAccordingToWatermarks
{
public:
    typedef flag_Type (*flagPolicy_ptr) ( flag_Type const&, flag_Type const& );

    //! The constructor
    /*!
     * @param flagToSet the flag we want to set if the condition is satisfied
     * @param watermarks The std::vector of watermarks to be used in the comparison
     * @param flagPolicy The policy we want to use. Default is turnOn, which
     * turns on the bit-flags defined in flagToSet. Other possibilities are turnOff or
     * change (the names are self explanatory)
     */
    SetFlagAccordingToWatermarks ( const flag_Type& flagToSet,
                                   const std::vector<markerID_Type>& watermarks,
                                   const flagPolicy_ptr& flagPolicy = &Flag::turnOn );

    //! The operator doing the job
    /*!
     * It is the operator that does the check on a given mesh entity
     * @param e The entity to change (possibly) using the given policy
     */
    template<typename MeshEntity>
    void operator() ( MeshEntity& e ) const
    {
        if ( std::binary_search ( M_watermarks.begin(), M_watermarks.end(), e.markerID() ) )
        {
            e.replaceFlag ( M_flagPolicy ( e.flag(), M_flagToSet ) );
        }
    }

private:
    const flag_Type M_flagToSet;
    std::vector<markerID_Type> M_watermarks;
    const flagPolicy_ptr M_flagPolicy;
};

/** @}*/
/** @defgroup markerIDchangers
 * Utility to change the marker ids according to some policy
 * @{
 */
//! Set markers according to a map.
/**
 * This utility is used in some LifeV solvers to change the marker IDs of certain entities on the
 * fly
 *
 * @param locDof Contains a map of integers whose second entry contains the entity id to be
 * changed
 * @param newMarker the new marker id
 *
 * @note  it was originally a method of Regionmesh3D named edgeMarkers. It has been generalised
 * and taken away from RegionMesh
 *
 * @todo Take away from here and put it in another header file (MeshUtility.hpp for instance)
 */
template <typename MeshEntityList>
void
ChangeMarkersAccordingToMap ( MeshEntityList& entityList,
                              std::map<UInt, UInt> const& locDof,
                              UInt const newMarker )
{
    typedef std::map<UInt, UInt>::const_iterator it_type;
    for ( it_type it = locDof.begin(); it != locDof.end(); ++it )
    {
        entityList[ it->second ].setMarkerID ( newMarker );
    }
}
/** @}*/

} // Namespace LifeV

#endif /* SELECTMARKER_H */
