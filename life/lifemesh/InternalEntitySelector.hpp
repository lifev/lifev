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
#include <life/lifemesh/Marker.hpp>
#include <life/lifearray/MeshEntityContainer.hpp>

namespace LifeV
{

/** InternalEntitySelector.
 *
 *  Functor class that tells whether an entity flag corresponds to an internal face
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
    InternalEntitySelector(const markerID_Type & w);
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
    bool operator()(markerID_Type const & test) const;
    //@}

private:

    //! The current watermark
    markerID_Type M_watermarkFlag;

}; // Class InternalEntitySelector


/** @defgroup MeshEntityFlagsFunctors
 * Useful functors to change entityflags according to certain conditions
 * @{
 */
/** Set the flag of a mesh entity according to the value of its Marker id.
 *
 *  We compare the marker id with a watermark according to a policy which is a comparison operator
 *  with signature bool operator()(markerID_Type const & mId, markerID_Type const & watermark)
 *  which by default is greater<>  (i.e. returns true if  mId > watermark).
 *  If the comparison is true the given flag is assigned to the entity using a second policy
 *  implemented as a function:
 *
 *  flag_Type flagPolicy  ( flag_Type const & inputFlag, flag_Type const & refFlag )
 *
 *  passed in the constructor and defaulted to Flag::turnOn
 *
 *  Example:
 *  I want to set the flag INTERNAL_INTERFACE to all faces with entity flag > 1000
 *
 * @verbatim
    SetFlagAccordingToWatermark<face_Type>
           changer(EntityFlags::INTERNAL_INTERFACE,1000,Flag::turnOn)
    // the last argument (turnOn) is not needed
    mesh.faceList.changeAccordingToFunctor(changer);
   @endverbatim
 *  I want to set the flag INTERNAL_INTERFACE to all faces with entity flag =12000
 *  @verbatim
 *  SetFlagAccordingToWatermark<face_Type,std::equal_to<markerID_Type> >
 *                          changer(EntityFlags::INTERNAL_INTERFACE,12000,Flag::turnOn)
 *  (the last argument is not needed)
 *  mesh.faceList.changeAccordingToFunctor(changer);
 *  @endverbatim
 *
 */
template<typename MeshEntity, typename Policy=std::greater<markerID_Type> >
    class SetFlagAccordingToWatermark{
    public:
    typedef flag_Type (*flagPolicy_ptr)(flag_Type const &, flag_Type const &);
    SetFlagAccordingToWatermark(const flag_Type & flagToSet, const markerID_Type & watermark,
                                const flagPolicy_ptr & flagPolicy=&Flag::turnOn ):
        M_flagToSet(flagToSet),M_watermark(watermark), M_policy(Policy()),M_flagPolicy(flagPolicy)
    {
    }
    void operator()(MeshEntity & e)const
    {
        if ( M_policy(e.marker(),M_watermark) ) e.setFlag(M_flagPolicy(e.flag(),M_flagToSet));
    }
    private:
    const flag_Type M_flagToSet;
    const markerID_Type M_watermark;
    const Policy M_policy;
    const flagPolicy_ptr M_flagPolicy;
};

/** Set the boolean flag of a mesh entity according to the value of a vector of Marker ids.
 *
 *  We compare the marker ids with all values contained in the vector
 *  using a policy which is a comparison operator
 *  with signature bool operator()(markerID_Type const & mId, markerID_Type const & watermark)
 *  which by default is equal_to<>  (i.e. returns true if  mId == watermark).
 *  If the comparison is true the given flag is assigned to the entity using a second policy implemented as a function
 *  @verbatim
 *  flag_Type flagPolicy  ( flag_Type const & inputFlag, flag_Type const & refFlag )
 *  @endverbatim
 *  passed in the constructor and defaulted to Flag::turnOn
 *
 *  Example:
 *  I want to set the flag INTERNAL_INTERFACE to all faces with entity flag = 1000, 2000 and 3000
 *  @verbatim
   vector<markerID_Type> fl; fl.push_back(1000); fl.push_back(2000); fl.push_back(3000)
   SetFlagAccordingToWatermarks<face_Type> changer(INTERNAL_INTERFACE,fl,Flag::turnOn)
   //the last argument is not needed
   mesh.faceList.changeAccordingToFunctor(changer);
   @endverbatim
 *
 */

template<typename MeshEntity, typename Policy=std::equal_to<markerID_Type> >
    class SetFlagAccordingToWatermarks{
    public:
    typedef flag_Type (*flagPolicy_ptr)(flag_Type const &, flag_Type const &);
    SetFlagAccordingToWatermarks(const flag_Type & flagToSet,
                                     const std::vector<markerID_Type> & watermarks,
                                     const flagPolicy_ptr & flagPolicy=&Flag::turnOn ):
        M_flagToSet(flagToSet),
        M_watermarks(watermarks),
        M_policy(Policy()),
        M_flagPolicy(flagPolicy){}
    void operator()(MeshEntity & e)const
    {
        typedef std::vector<markerID_Type>::iterator it;
        for (it i=M_watermarks.begin();i<M_watermarks.end();++i)
            if ( M_policy( e.marker(),*i ) ){
                e.setFlag( M_flagPolicy(e.flag(),M_flagToSet) );
                break;
            }
    }
    private:
    const flag_Type M_flagToSet;
    const std::vector<markerID_Type> M_watermarks;
    const Policy M_policy;
    const flagPolicy_ptr M_flagPolicy;
};

/** @}*/
/** @defgroup markerIDchangers
 * Utility to change the marker ids according to some policy
 * @{
 */
//! Set markers according to a map.
/**
 * This utility is used in some LifeV solvers to change the markers of certain entities on the
 * fly
 *
 * @param locDof Contains a map of int pairs whose second entry contains the entity id to be
 * changed
 * @param TimeAdvanceNewMarker the new marker id
 *
 * @note, it was originally a method of regionemesh3D names edgeMarkers. It has been generalised
 * and taken away from regionmesh
 */
template <typename MeshEntityList>
void
ChangeMarkersAccordingToMap(MeshEntityList & entityList,
                            std::map<UInt, UInt> const& locDof,
                            UInt const TimeAdvanceNewmarker)
{
    typedef std::map<UInt, UInt>::const_iterator it_type;
    for (it_type IT=locDof.begin(); IT!=locDof.end(); ++IT)
        entityList[IT->second].setMarker(TimeAdvanceNewmarker);
}
/** @}*/
} // Namespace LifeV

#endif /* SELECTMARKER_H */
