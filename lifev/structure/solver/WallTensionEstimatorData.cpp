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
 *  @file
 *  @brief DataStructure - File containing a data container for wall tension analysis
 *
 *  @version 1.0
 *  @date 19-04-2012
 *  @author Paolo Tricerri
 *
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#include <lifev/structure/solver/WallTensionEstimatorData.hpp>

namespace LifeV
{

//=====================================================
// Constructors
//=====================================================
WallTensionEstimatorData::WallTensionEstimatorData() :
    M_nameFile                         ( ),
    M_analysisType                     ( ),
    M_recoveryVariable                 ( ),
    M_initialTime                      ( ),
    M_finalTime                        ( ),
    M_iterStart                        ( ),
    M_iterEnd                          ( )
{
}

WallTensionEstimatorData::WallTensionEstimatorData ( const WallTensionEstimatorData& wallTensionEstimatorData ) :
    M_nameFile                         ( wallTensionEstimatorData.M_nameFile ),
    M_analysisType                     ( wallTensionEstimatorData.M_analysisType ),
    M_recoveryVariable                 ( wallTensionEstimatorData.M_recoveryVariable ),
    M_initialTime                      ( wallTensionEstimatorData.M_initialTime ),
    M_finalTime                        ( wallTensionEstimatorData.M_finalTime ),
    M_iterStart                        ( wallTensionEstimatorData.M_iterStart ),
    M_iterEnd                          ( wallTensionEstimatorData.M_iterEnd )
{
}

// ===================================================
// Operators
// ===================================================
WallTensionEstimatorData&
WallTensionEstimatorData::operator= ( const WallTensionEstimatorData& wallTensionEstimatorData )
{
    if ( this != &wallTensionEstimatorData )
    {
        M_nameFile                         = wallTensionEstimatorData.M_nameFile;
        M_analysisType                     = wallTensionEstimatorData.M_analysisType;
        M_recoveryVariable                 = wallTensionEstimatorData.M_recoveryVariable;
        M_initialTime                      = wallTensionEstimatorData.M_initialTime;
        M_finalTime                        = wallTensionEstimatorData.M_finalTime;
        M_iterStart                        = wallTensionEstimatorData.M_iterStart;
        M_iterEnd                          = wallTensionEstimatorData.M_iterEnd;

    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
WallTensionEstimatorData::setup ( const GetPot& dataFile, const std::string& section )
{

    // physics
    M_nameFile = dataFile ( ( section + "/analysis/nameFile" ).data(), "NO_DEFAULT_NAME_FILE" );
    M_typeFile = dataFile ( ( section + "/analysis/typeFile" ).data(), "NO_DEFAULT_FILE_TYPE" );

    M_analysisType = dataFile ( ( section + "/analysis/analysisType" ).data(), "NO_DEFAULT_ANALYSIS_TYPE" );
    M_recoveryVariable = dataFile ( ( section + "/analysis/recoveryVariable" ).data(), "NO_DEFAULT_ANALYSIS_TYPE" );

    UInt timesNumber (0);

    if ( !M_analysisType.compare ("istant") )
    {
        timesNumber = 1;
    }
    else
    {
        timesNumber = dataFile.vector_variable_size ( ( section + "/analysis/start" ).data() );
    }

    //Resizing the vectors to store the right number of iterations & times
    M_initialTime.resize (timesNumber);
    M_finalTime.resize (timesNumber);

    M_iterStart.resize (timesNumber);
    M_iterEnd.resize (timesNumber);

    for ( UInt i (0); i < timesNumber ; ++i)
    {
        M_initialTime[i] = dataFile ( ( section + "/analysis/start" ).data(), 0., i );
        M_finalTime[i]   = dataFile ( ( section + "/analysis/end"   ).data(), 0., i );

        M_iterStart[i] = dataFile ( ( section + "/analysis/iterationStart" ).data(), "00000", i );
        M_iterEnd[i]   = dataFile ( ( section + "/analysis/iterationEnd"   ).data(), "00000", i );
    }
}

void
WallTensionEstimatorData::showMe ( std::ostream& output ) const
{
    // physics
    output << "\n*** Values for wall tension analysis [solid/analysis]\n\n";
    output << "Name File                = " << M_nameFile << std::endl;
    output << "Type File                = " << M_typeFile << std::endl;
    output << "Analysis Type            = " << M_analysisType << std::endl;
    output << "Recovery Variable        = " << M_recoveryVariable << std::endl;
    output << "The numbers of intervals is =  " << M_initialTime.size() << std::endl;

    for ( UInt i (0); i < M_initialTime.size() ; ++i )
    {
        output <<  i + 1 << "  Interval: " <<  std::endl;
        output << "initial Time " << i + 1 << "  = " << M_initialTime[i] << std::endl;
        output << "final  Time " << i + 1 << "   = " << M_finalTime[i] << std::endl;

        output << "iteration Start " << i + 1 << "= " << M_iterStart[i] << std::endl;
        output << "iteration End " << i + 1 << "= " << M_iterEnd[i] << std::endl;
    }
}

}
