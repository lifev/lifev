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

#include <lifev/electrophysiology/stimulus/StimulusPacingProtocol.hpp>

namespace LifeV
{

// ===================================================
//! Constructors
// ===================================================
StimulusPacingProtocol::StimulusPacingProtocol() :
    M_radius ( 0 ),
    M_stimulusAmplitude ( 0 ),
    M_pacingSite_X ( 0 ),
    M_pacingSite_Y ( 0 ),
    M_pacingSite_Z ( 0 ),
    M_numberStimulus ( 0 ),
    M_nbStimMax(-1)
{
}

// ===================================================
//! Methods
// ===================================================


Real StimulusPacingProtocol::pacingProtocolChoice ( const Real& t)
{

    if ( M_pacingProtocol == "FCL" )
    {
        return fixedCycleLength ( t );
    }

    else if ( M_pacingProtocol == "FCL-ExtraSt" )
    {
        return fixedCycleLengthwExtraStim ( t );
    }

    else if ( M_pacingProtocol == "S1S2" )
    {
        return standardS1S2Protocol ( t );
    }

    else if ( M_pacingProtocol == "DynPro" )
    {
        return dynamicProtocol ( t );
    }

    else
    {
        return fixedCycleLength ( t );
    }


}


Real StimulusPacingProtocol::fixedCycleLength ( const Real& t )
{
    Real current = 0;
    if ( M_numberStimulus < M_nbStimMax || M_nbStimMax == -1 )
    {
        if ( t >= M_startingTimeStimulus && t <= M_startingTimeStimulus + M_stimulusDuration )
        {
            current = M_stimulusAmplitude;
            if ( t + M_dt >=  M_startingTimeStimulus + M_stimulusDuration )
            {
            	std::cout << "\nUpdating the stimulus interval";
                M_numberStimulus++;
                M_startingTimeStimulus = M_startingTimeStimulus + M_stimulusInterval;
            }
        }
        else
        {
            current = 0;
        }
    }
    else
    {
        current = 0;
    }
    return current;
}

Real StimulusPacingProtocol::fixedCycleLengthwExtraStim ( const Real& t )
{
    Real current = 0;
    if ( M_numberStimulus < M_nbStimMax || M_nbStimMax == -1 )
    {
        if ( t >= M_startingTimeStimulus && t <= M_startingTimeStimulus + M_stimulusDuration )
        {
            current = M_stimulusAmplitude;

            if ( t + M_dt >=  M_startingTimeStimulus + M_stimulusDuration  )
            {
                if ( M_numberStimulus < M_repeatSt )
                {
                    M_startingTimeStimulus = M_startingTimeStimulus + M_stimulusInterval;
                    M_numberStimulus++;
                }
                else
                {
                    if ( M_pacingProtocolType == "S1-S2" )
                    {
                        M_startingTimeStimulus = M_startingTimeStimulus + M_stIntS1S2;
                        M_numberStimulus = 0;
                    }
                    else if ( M_pacingProtocolType == "S1-S2-S3" )
                    {
                        if ( M_numberStimulus == M_repeatSt )
                        {
                            M_startingTimeStimulus = M_startingTimeStimulus + M_stIntS1S2;
                            M_numberStimulus++;
                        }
                        else
                        {
                            M_startingTimeStimulus = M_startingTimeStimulus + M_stIntS2S3;
                            M_numberStimulus = 0;
                        }
                    }
                    else if ( M_pacingProtocolType == "S1-S2-S3-S4" )
                    {
                        if ( M_numberStimulus == M_repeatSt )
                        {
                            M_startingTimeStimulus = M_startingTimeStimulus + M_stIntS1S2;
                            M_numberStimulus++;
                        }
                        else if ( M_numberStimulus == M_repeatSt + 1 )
                        {
                            M_startingTimeStimulus = M_startingTimeStimulus + M_stIntS2S3;
                            M_numberStimulus++;
                        }
                        else
                        {
                            M_startingTimeStimulus = M_startingTimeStimulus + M_stIntS3S4;
                            M_numberStimulus = 0;
                        }
                    }
                    else
                    {
                        M_startingTimeStimulus = M_startingTimeStimulus + M_stimulusInterval;
                    }
                }
            }
        }
    }

    return current;
}

Real StimulusPacingProtocol::standardS1S2Protocol ( const Real& t)
{
    Real current = 0;
    if ( t < M_nbStimMax * M_stimulusInterval )
    {
        if ( t >= M_startingTimeStimulus && t <= M_startingTimeStimulus + M_stimulusDuration )
        {
            current = M_stimulusAmplitude;

            if ( t + M_dt >=  M_startingTimeStimulus + M_stimulusDuration )
            {
                M_startingTimeStimulus = M_startingTimeStimulus + M_stimulusInterval;

                if ( t > ( M_nbStimMax - 1 ) * M_stimulusInterval && t < M_nbStimMax * M_stimulusInterval )
                {
                    M_numberStimulus = 0;
                }
                else
                {
                    M_numberStimulus++;
                }
            }
        }
    }
    else
    {
        if ( M_stIntS1S2 >= M_stIntS1S2Min )
        {
            if ( t >= M_startingTimeStimulus && t <= M_startingTimeStimulus + M_stimulusDuration)
            {
                current = M_stimulusAmplitude;

                if ( t >= M_startingTimeStimulus + M_stimulusDuration - M_dt && t <= M_startingTimeStimulus + M_stimulusDuration )
                {
                    M_numberStimulus++;

                    if ( M_numberStimulus < M_repeatSt )
                    {
                        M_startingTimeStimulus = M_startingTimeStimulus + M_stimulusInterval;
                    }

                    else if ( M_numberStimulus == M_repeatSt )
                    {
                        M_startingTimeStimulus = M_startingTimeStimulus + M_stIntS1S2;
                    }

                    else if ( M_numberStimulus == M_repeatSt + 1 )
                    {
                        M_startingTimeStimulus = M_startingTimeStimulus + M_stimulusInterval;
                        M_numberStimulus = 0;

                        if ( M_stIntS1S2 > 100 )
                        {
                            M_stIntS1S2 = M_stIntS1S2 - 100;
                        }

                        else if ( M_stIntS1S2 <= 100 && M_stIntS1S2 > 30 )
                        {
                            M_stIntS1S2 = M_stIntS1S2 - 5;
                        }

                        else if (  M_stIntS1S2 <= 30 &&  M_stIntS1S2 > 20 )
                        {
                            M_stIntS1S2 = M_stIntS1S2 - 1;
                        }

                        else if ( M_stIntS1S2 <= 20 )
                        {
                            M_stIntS1S2 = M_stIntS1S2 - 5;
                        }
                    }
                }
            }
        }
    }
    return current;
}

Real StimulusPacingProtocol::dynamicProtocol ( const Real& t )
{
    Real current = 0;
    if ( M_stimulusInterval >= M_minimumStimulusInterval )
    {
        if ( t >= M_startingTimeStimulus && t <= M_startingTimeStimulus + M_stimulusDuration )
        {
            current = M_stimulusAmplitude;

            if ( t >= M_startingTimeStimulus + M_stimulusDuration - M_dt && t <= M_startingTimeStimulus + M_stimulusDuration )
            {
                M_numberStimulus++;
                M_startingTimeStimulus = M_startingTimeStimulus + M_stimulusInterval;
            }
        }
        else
        {
            current = 0;
        }

        if ( t > M_tShortS1S1 )
        {
            if ( M_stimulusInterval > 1000 )
            {
                M_stimulusInterval      = M_stimulusInterval - 1000;
                M_tShortS1S1 = M_tShortS1S1 + M_stimulusInterval * 20;
            }
            else if ( M_stimulusInterval <= 1000 && M_stimulusInterval > 300 )
            {
                M_stimulusInterval      = M_stimulusInterval - 50;
                M_tShortS1S1 = M_tShortS1S1 + M_stimulusInterval * 20;
            }
            else if (  M_stimulusInterval <= 300 &&  M_stimulusInterval > 200 )
            {
                M_stimulusInterval      = M_stimulusInterval - 10;
                M_tShortS1S1 = M_tShortS1S1 + M_stimulusInterval * 20;
            }
            else if ( M_stimulusInterval <= 200 )
            {
                M_stimulusInterval      = M_stimulusInterval - 5;
                M_tShortS1S1 = M_tShortS1S1 + M_stimulusInterval * 20;
            }
        }
    }
    else
    {
        current = 0;
    }
    return current;
}

Real StimulusPacingProtocol::appliedCurrent ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*i*/ )
{

    Real current = 0.0;
    const Real volumeOfBall = (4. / 3.) * M_PI * M_radius * M_radius * M_radius;
    Real distance = std::sqrt ( (x - M_pacingSite_X) * (x - M_pacingSite_X) + (y - M_pacingSite_Y) * (y - M_pacingSite_Y) + (z - M_pacingSite_Z) * (z - M_pacingSite_Z) );
    if (distance <= M_radius )
    {
       current += pacingProtocolChoice( t ) / volumeOfBall;
    }
    return current;

}

void StimulusPacingProtocol::showMe()
{
    std::cout << "\n\n\t\tPacing protocol Informations\n\n";

    std::cout << "\n\t\tList of parameters:\n\n";

    std::cout << "Istim: " << M_stimulusAmplitude << std::endl;
    std::cout << "StimDuration: " << M_stimulusDuration << std::endl;
    std::cout << "1st stimuli time: " << M_startingTimeStimulus << std::endl;
    std::cout << "Pacing protocol: " << M_pacingProtocol << std::endl;
    std::cout << "Pacing site: " << M_pacingSite_X << " " << M_pacingSite_Y << " " << M_pacingSite_Z << std::endl;
    std::cout << "Radius stimulus: " << M_radius << std::endl;
    if ( M_pacingProtocol == "FCL" )
    {
        std::cout << "S1-S1 interval: " << M_stimulusInterval << std::endl;
        std::cout << "NbStimuliMax: " << M_nbStimMax << std::endl;
    }
    else if ( M_pacingProtocol == "FCL-ExtraSt" )
    {
        std::cout << "Pacing protocol type: " << M_pacingProtocolType << std::endl;
        std::cout << "S1-S1 interval: " << M_stimulusInterval << std::endl;
        std::cout << "NbStimuliMax: " << M_nbStimMax << std::endl;
        std::cout << "Repeat S1 stimuli: " << M_repeatSt << std::endl;

        if ( M_pacingProtocolType == "S1-S2" )
        {
            std::cout << "S1-S2 interval: " << M_stIntS1S2 << std::endl;
        }
        else if ( M_pacingProtocolType == "S1-S2-S3" )
        {
            std::cout << "S1-S2 interval: " << M_stIntS1S2 << std::endl;
            std::cout << "S2-S3 interval: " << M_stIntS2S3 << std::endl;
        }
        else if ( M_pacingProtocolType == "S1-S2-S3-S4" )
        {
            std::cout << "S1-S2 interval: " << M_stIntS1S2 << std::endl;
            std::cout << "S2-S3 interval: " << M_stIntS2S3 << std::endl;
            std::cout << "S3-S4 interval: " << M_stIntS3S4 << std::endl;
        }

    }
    else if ( M_pacingProtocol == "S1S2Pro" )
    {
        std::cout << "S1-S1 interval: " << M_stimulusInterval << std::endl;
        std::cout << "NbStimuliMax for stabilisation: " << M_nbStimMax << std::endl;
        std::cout << "First S1-S2 interval: " << M_stIntS1S2 << std::endl;
        std::cout << "Minimum S1-S2 interval: " << M_stIntS1S2Min << std::endl;
        std::cout << "Repeat S1 stimuli: " << M_repeatSt << std::endl;
    }
    else if ( M_pacingProtocol == "DynPro" )
    {
        std::cout << "First S1-S1 interval: " << M_stimulusInterval << std::endl;
        std::cout << "Minimum S1-S1 interval: " << M_minimumStimulusInterval << std::endl;
        std::cout << "First time S1-S1 interval decrease: " << M_tShortS1S1 << std::endl;
    }
    else
    {
        std::cout << "S1-S1 interval: " << M_stimulusInterval << std::endl;
        std::cout << "NbStimuliMax: " << M_nbStimMax << std::endl;
    }
    std::cout << "\n\t\t End of Pacing protocol Informations\n\n\n";
}

}
