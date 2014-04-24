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
 @brief Class for applying cardiac stimulus represented by
	a current at a single point with a given time and duration

 @date 02-2014
 @author Simone Palamara <palamara.simone@gmail.com>

 @last update 02-2014
 */

#ifndef STIMULUSSINGLESOURCE_HPP_
#define STIMULUSSINGLESOURCE_HPP_

#include <lifev/electrophysiology/stimulus/ElectroStimulus.hpp>

namespace LifeV
{

class StimulusSingleSource : public ElectroStimulus
{

public:


    //! @name Constructors & Destructor
    //@{

    //!Empty Constructor
    /*!
     */
    StimulusSingleSource();

    //! Destructor
    virtual ~StimulusSingleSource() {}

    //@}

    //! @name Set Methods
    //@{

    inline void setRadius ( Real r )
    {
        ASSERT (r > 0, "Invalid radius value.");
        M_radius = r;
    }

    inline void setTotalCurrent ( Real I )
    {
        ASSERT (I >= 0, "Invalid current value.");
        M_totalCurrent = I;
    }

    inline void setPacingSite ( Real x , Real y , Real z )
    {
        M_pacingSite_X = x;
        M_pacingSite_Y = y;
        M_pacingSite_Z = z;
    }

    inline void setStimDuration ( Real duration )
    {
        ASSERT (duration>0, "Invalid stimulus duration value.");
        M_StimDuration=duration;
    }

    inline void setStartingTimeStimulus ( Real startingTimeStimulus )
    {
        ASSERT (startingTimeStimulus>=0, "Invalid starting time stimulus.");
        M_startingTimeStimulus=startingTimeStimulus;
    }


	void setParameters (list_Type&  list)
	{
		this->setRadius ( list.get ( "applied_current_radius", 0.2 ) );
		this->setTotalCurrent ( list.get ( "applied_total_current", 1.0 ) );
        this->setPacingSite (list.get ( "pacing_site_X", 1.0 ),list.get ( "pacing_site_Y", 1.0 ),list.get ( "pacing_site_Z", 1.0 ) );
        this->setStimDuration( list.get ( "duration_stimulus", 1.0 ) );
        this->setStartingTimeStimulus( list.get ("starting_time_stimulus", 0.0 ) );
    }


    //@}

    //! @name Copy Methods
    //@{

    //@}

    //! @name Methods
    //@{
    Real appliedCurrent ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i );

    void showMe ();
    //@}

private:

    Real         M_radius;
    Real         M_totalCurrent;
    Real 		 M_pacingSite_X;
    Real 		 M_pacingSite_Y;
    Real  		 M_pacingSite_Z;
    Real 		 M_startingTimeStimulus;
    Real   		 M_StimDuration;
};

} // namespace LifeV

#endif /* STIMULUSSINGLESOURCE_HPP_ */
