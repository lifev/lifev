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

#ifndef WallTensionEstimatorData_H
#define WallTensionEstimatorData_H

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

// STL classes
#include <string>
#include <iostream>
#include <map>

// Boost classes
#include <boost/shared_ptr.hpp>

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LifeV core includes
#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/StringUtility.hpp>
#include <lifev/core/filter/GetPot.hpp>

namespace LifeV
{

//! DataElasticStructure - Data container for solid problems with elastic structure
class WallTensionEstimatorData
{
public:

    //! @name Type definitions
    //@{

    typedef  std::vector<std::string > iteration_Type;
    typedef  std::vector< Real >       time_Type;

    //@}


    //! @name Constructors & Destructor
    //@{0

    //! Empty Constructor
    WallTensionEstimatorData();

    //! Copy constructor
    /*!
     * @param WallTensionEstimatorData - WallTensionEstimatorData
     */
    WallTensionEstimatorData ( const WallTensionEstimatorData& wallTensionEstimatorData );

    //! Destructor
    virtual ~WallTensionEstimatorData() {};

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param WallTensionEstimatorData - WallTensionEstimatorData
     */
    WallTensionEstimatorData& operator= ( const WallTensionEstimatorData& wallTensionEstimatorData );

    //@}


    //! @name Methods
    //@{

    //! Read the dataFile and set all the quantities
    /*!
     * @param dataFile data file
     * @param section section of the file
     */
    void setup ( const GetPot& dataFile, const std::string& section = "solid" );

    //! Display the values
    void showMe ( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set methods
    //@{

    //! Set nameFile
    /*!
     * @param nameFile name of the post-processing file to be opened
     */
    void setNameFile ( const  std::string& nameFile)
    {
        M_nameFile = nameFile;
    }

    //! Set Analysis Type
    /*!
     * @param analysis Type of analysis that can be performed
     */
    void setAnalysisType ( const std::string& analysisType)
    {
        M_analysisType = analysisType;
    }

    //! Set Recovery Variable
    /*!
     * @param recoveryVariable the field that is recovered (displacement or tensions)
     */
    void setRecoveryVariable ( const std::string& recoveryVariable)
    {
        M_recoveryVariable = recoveryVariable;
    }

    //! Set initialTime
    /*!
     * @param initial Time initial of the analysis
     */
    void setInitialTime (const UInt i, const Real initialTime )
    {
        M_initialTime[i] = initialTime;
    }

    //! Set final Time
    /*!
     * @param final Time final time of the analysis
     */
    void setFinalTime ( const UInt i, const Real finalTime )
    {
        M_finalTime[i] = finalTime;
    }

    //! Set starting iteration
    /*!
     * @param starting iteration
     */
    void setIterationStart ( const UInt i, const std::string& iterStart )
    {
        M_iterStart[i] = iterStart;
    }

    //! Set final Time
    /*!
     * @param final Time final time of the analysis
     */
    void setIterationEnd ( const UInt i, const std::string& iterEnd )
    {
        M_iterEnd[i] = iterEnd;
    }

    //@}


    //! @name Get methods
    //@{

    //! Get nameFile
    /*!
     * @return std::string with the name of the file used for the analysis
     */
    const std::string& nameFile() const
    {
        return M_nameFile;
    }

    //! Get typeFile
    /*!
     * @return std::string with the type of the file used for the analysis
     */
    const std::string& typeFile() const
    {
        return M_typeFile;
    }

    //! Get analysisType
    /*!
     * @return std::string with the type of analysis that has to be performed
     */
    const std::string& analysisType() const
    {
        return M_analysisType;
    }

    //! Get recoveryVariable
    /*!
     * @return std::string with the name of the field that is recovered
     */
    const std::string& recoveryVariable() const
    {
        return M_recoveryVariable;
    }

    //! Get initial Time
    /*!
     * @return initial Time of the analysis
     */
    const time_Type& initialTime() const
    {
        return M_initialTime;
    }

    //! Get initial Time of the interval i
    /*!
     * @return initial Time of the i-th interval
     */
    const Real initialTime (const UInt i) const
    {
        return M_initialTime[i];
    }

    //! Get final Time
    /*!
     * @return final time of the anaysis
     */
    const time_Type& finalTime() const
    {
        return M_finalTime;
    }

    //! Get final Time of the i-th interval
    /*!
     * @return final time of the i-th interval
     */
    const Real finalTime (const UInt i) const
    {
        return M_finalTime[i];
    }

    //! Get starting iteration
    /*!
     * @return starting iteration
     */
    const iteration_Type& iterStart() const
    {
        return M_iterStart;
    }

    //! Get starting iteration of the i-th interval
    /*!
     * @return starting iteration of the i-th interval
     */
    const std::string& iterStart (const UInt i) const
    {
        return M_iterStart[i];
    }


    //! Get ending iteration
    /*!
     * @return ending iteration
     */
    const iteration_Type& iterEnd() const
    {
        return M_iterEnd;
    }

    //! Get starting iteration
    /*!
     * @return starting iteration
     */
    const std::string& iterEnd (const UInt i) const
    {
        return M_iterEnd[i];
    }

    //@}

private:

    enum M_analysisType {istant, interval};
    //! Name of the file
    std::string           M_nameFile;
    std::string           M_typeFile;

    //! Type of Analysis
    std::string           M_analysisType;
    //! Field that is recovered
    std::string           M_recoveryVariable;

    //! Initial & Final Time
    time_Type             M_initialTime;
    time_Type             M_finalTime;

    //! Initial & Final iteration
    iteration_Type        M_iterStart;
    iteration_Type        M_iterEnd;


};

} // end namespace LifeV

#endif // WallTensionEstimatorData_H
