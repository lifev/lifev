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
 @brief Class for applying cardiac stimulus at a single point according to a pacing protocol

 @date 02-2014
 @author Simone Palamara <palamara.simone@gmail.com>

 @last update 02-2014
 */

#ifndef STIMULUSPACINGPROTOCOL_HPP_
#define STIMULUSPACINGPROTOCOL_HPP_

#include <lifev/electrophysiology/stimulus/ElectroStimulus.hpp>

namespace LifeV
{

class StimulusPacingProtocol : public ElectroStimulus
{

public:

    //! @name Constructors & Destructor
    //@{

    //!Empty Constructor
    /*!
     */
    StimulusPacingProtocol();

    //! Destructor
    virtual ~StimulusPacingProtocol() {}

    //@}

    //! @name Set Methods
    //@{

    inline void setRadius ( Real r )
    {
        ASSERT (r > 0, "Invalid radius value.");
        M_radius = r;
    }

    inline void setStimulusAmplitude ( Real I )
    {
        ASSERT (I >= 0, "Invalid current value.");
        M_stimulusAmplitude = I;
    }

    inline void setPacingSite ( Real x , Real y , Real z )
    {
        M_pacingSite_X = x;
        M_pacingSite_Y = y;
        M_pacingSite_Z = z;
    }

    inline void setPacingProtocol ( std::string PacingProtocol )
    {
        M_pacingProtocol = PacingProtocol;
    }


    inline const Real& StInt() const
    {
        return M_stimulusInterval;
    }
    inline void setStInt (const Real& StInt)
    {
        this->M_stimulusInterval = StInt;
    }

    inline const Real& timeShortS1S1() const
    {
        return M_tShortS1S1;
    }
    inline void setTimeShortS1S1 (const Real& tShortS1S1)
    {
        this->M_tShortS1S1 = tShortS1S1;
    }

    inline const Real& startingTimeStimulus() const
    {
        return M_startingTimeStimulus;
    }
    inline void setStartingTimeStimulus (const Real& startingTimeStimulus)
    {
        this->M_startingTimeStimulus = startingTimeStimulus;
    }

    inline const Real& stIntMin() const
    {
        return M_minimumStimulusInterval;
    }
    inline void setStIntMin (const Real& stIntMin)
    {
        this->M_minimumStimulusInterval = stIntMin;
    }

    inline const Real& timeStep() const
    {
        return M_dt;
    }
    inline void setTimeStep (const Real& dt)
    {
        this->M_dt = dt;
    }


    inline const Real& stIntS1S2() const
    {
        return M_stIntS1S2;
    }
    inline void setStIntS1S2 (const Real& stIntS1S2 )
    {
        this->M_stIntS1S2 = stIntS1S2;
    }

    inline const Real& stIntS1S2Min() const
    {
        return M_stIntS1S2Min;
    }
    inline void setStIntS1S2Min (const Real& stIntS1S2Min )
    {
        this->M_stIntS1S2Min = stIntS1S2Min;
    }

    inline const Real& stIntS2S3() const
    {
        return M_stIntS2S3;
    }
    inline void setStIntS2S3 (const Real& stIntS2S3)
    {
        this->M_stIntS2S3 = stIntS2S3;
    }

    inline const Real& stIntS3S4() const
    {
        return M_stIntS3S4;
    }
    inline void setStIntS3S4 (const Real& stIntS3S4)
    {
        this->M_stIntS3S4 = stIntS3S4;
    }

    inline const int& nbStimMax() const
    {
        return M_nbStimMax;
    }
    inline void setNbStimMax (const int& nbStimMax)
    {
        this->M_nbStimMax = nbStimMax;
    }

    inline const int& repeatSt() const
    {
        return M_repeatSt;
    }
    inline void setRepeatSt (const int& repeatSt)
    {
        this->M_repeatSt = repeatSt;
    }

    inline const Real& stimulusDuration() const
    {
        return M_stimulusDuration;
    }
    inline void setStimulusDuration (const Real& stimDur)
    {
        this->M_stimulusDuration = stimDur;
    }

    inline void setPacingProtocolType (const std::string& pacProTyp)
    {
        this->M_pacingProtocolType = pacProTyp;
    }


    void setParameters (list_Type&  list)
    {
        M_startingTimeStimulus =  list.get ("starting_time_stimulus", 0.0);
        M_tShortS1S1 = list.get ("tShortS1S1", 0.001);
        M_stimulusInterval = list.get ("stimulus_interval", 0.0 );
        this->setStIntMin ( list.get ("stIntMin", 0.001) );
        this->setStIntS1S2 ( list.get ("stIntS1S2", 0.0) );
        this->setStIntS1S2Min ( list.get ("stIntS1S2Min", 10.0) );
        this->setStIntS2S3 ( list.get ("stIntS2S3", 0.01) );
        this->setStIntS3S4 ( list.get ("stIntS3S4", 0.01) );
        M_nbStimMax = list.get ("nbStimMax", 1);
        this->setRepeatSt ( list.get ("repeatSt", 1) );
        this->setTimeStep ( list.get ("dt", 0.01) );
        this->setPacingProtocol ( list.get ("pacPro", "FCL-ExtraSt") );
        this->setPacingProtocolType ( list.get ("pacProType", "S1-S2") );
        M_stimulusDuration = list.get ("duration_stimulus", 1.0 );
        M_pacingSite_X = list.get ( "pacing_site_X", 1.0 );
        M_pacingSite_Y = list.get ( "pacing_site_Y", 1.0 );
        M_pacingSite_Z = list.get ( "pacing_site_Z", 1.0 );
        M_radius = list.get ( "applied_current_radius", 0.2 );
        M_stimulusAmplitude = list.get ( "stimulus_amplitude", 1.0 );
    }


    //@}

    //! @name Copy Methods
    //@{

    //@}

    //! @name Methods
    //@{
    Real appliedCurrent ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i );


    // Stimulation at constant frequency
    Real fixedCycleLength ( const Real& t );

    // Stimulation at constant frequency with extra stimulation ( respective Si-Sj interval = constant )
    Real fixedCycleLengthwExtraStim ( const Real& t );

    // Standard stimulation protocol
    Real standardS1S2Protocol ( const Real& t );

    // Dynamic stimulation protocol
    Real dynamicProtocol ( const Real& t );

    Real pacingProtocolChoice ( const Real& t );


    void showMe();

    //@}

private:

    Real         M_radius;
    Real         M_stimulusAmplitude;
    Real         M_pacingSite_X;
    Real         M_pacingSite_Y;
    Real         M_pacingSite_Z;

    // Values of the stimulation interval used in the protocols.
    Real M_startingTimeStimulus;
    Real M_tShortS1S1;
    Real M_stimulusInterval;                   //Duration of the time interval between two stimuli
    Real M_minimumStimulusInterval;                //Minimum duration of the time interval between two stimuli in a dynamic pacing protocol
    Real M_stIntS1S2;
    Real M_stIntS1S2Min;
    Real M_stIntS2S3;
    Real M_stIntS3S4;

    // Number of stimulation information
    int M_nbStimMax;                //Number of stimuli that we want to apply // infinity = -1
    int M_repeatSt;
    int M_numberStimulus;

    // Choice of the protocol
    std::string M_pacingProtocol;
    std::string M_pacingProtocolType;

    // Value of the current stimulation
    Real M_stimulusDuration;

    //Discretization time step used to understand when a stimulus in
    //the pacing train
    Real M_dt;



};

} // namespace LifeV

#endif /* STIMULUSPACINGPROTOCOL_HPP_ */
