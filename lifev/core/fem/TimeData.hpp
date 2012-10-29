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
    @brief File containing a class for handling temporal discretization

    @author M.A. Fernandez
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @date 01-06-2009

    @contributor Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>

    @maintainer Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
 */

#ifndef TimeData_H
#define TimeData_H 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <ostream>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/LifeV.hpp>

namespace LifeV
{

//! TimeData - Class for handling temporal discretization.
/*!
 *  @author Cristiano Malossi
 *
 *  The class is a data container for the time discretization.
 */
class TimeData
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    TimeData();

    //! Constructor
    /*!
     * @param dataFile GetPot data file
     * @param section the section on the data file that contains the information on the time discretization
     */
    TimeData( const GetPot& dataFile, const std::string& section = "time_discretization" );

    //! Copy constructor
    /*!
     * @param TimeData - TimeData class
     */
    TimeData( const TimeData& TimeData);

    //! Virtual destructor
    virtual ~TimeData() {}

    //@}


    //! @name Methods
    //@{

    //! Read the dataFile and set all the internal quantities
    /*!
     * @param dataFile data file
     * @param section section of the file
     */
    void setup( const GetPot& dfile, const std::string& section = "time_discretization" );

    //! Update the time by a timestep.
    void updateTime() { M_time += M_timeStep; ++M_timeStepNumber; }

    //! Return if we can make a new timestep
    /*!
     * @return true if time <= endTime, false viceversa
     */
    bool canAdvance() { return round( M_time ) <= round( M_endTime ); }

    //! Return if it is the initial time step
    /*!
     * @return true if time = initial time, false viceversa
     */
    bool isFirstTimeStep() { return round( M_time ) == round( M_initialTime ); }

    //! Return if it is the N time step with respect to the current initial time
    /*!
     * @param n time step number
     * @return true if time - (n - 1) * timeStep == initialTime, false viceversa
     */
    bool isFirstNTimeStep( const UInt& n ) { return round( M_time - (n - 1) * M_timeStep ) == round( M_initialTime ); }

    //! Return if it is the last time step
    /*!
     * @return true if time + timestep > endTime, false viceversa.
     */
    bool isLastTimeStep() { return round( M_time + M_timeStep ) > round( M_endTime ); }

    //! Display general information about the content of the class
    /*!
        @param output - specify the output format (std::cout by default)
     */
    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the initial time step
    /*!
     * @param initialTime initial time step value
     */
    void setInitialTime( const Real& initialTime ) { M_initialTime = initialTime; }

    //! Set the final time step
    /*!
     * @param endtime final time step value
     */
    void setEndTime( const Real& endTime ) { M_endTime = endTime; }

    //! Set the present time of the simulation
    /*!
     * @param time present time value
     */
    void setTime( const Real& time ) { M_time = time; }

    //! Set the initial time step
    /*!
     * @param timeStep initial time step value
     */
    void setTimeStep( const Real& timeStep ) { M_timeStep = timeStep; }

    //! Set the time step number
    /*!
     * @param timeStepNumber time step number
     */
    void setTimeStepNumber( const UInt& timeStepNumber ) { M_timeStepNumber = timeStepNumber; }

    //@}


    //! @name Get Methods
    //@{

    //! Get the initial time step
    /*!
     * @return initial time step value
     */
    const Real& initialTime() const { return M_initialTime; }

    //! Get the final time step
    /*!
     * @return final time step value
     */
    const Real& endTime() const { return M_endTime; }

    //! Get the period
    /*!
     * @return period value
     */
    const Real& periodTime() const { return M_periodTime; }

    //! Get the present time
    /*!
     * @return time value
     */
    const Real& time() const { return M_time; }

    //! Get the time left
    /*!
     * @return time left value
     */
    Real leftTime() const { return round( M_endTime - M_time ); }

    //! Get the elapsed time
    /*!
     * @return elapsed time value
     */
    Real elapsedTime() const { return round( M_time - M_initialTime ); }

    //! Get the present time shifted inside the first cycle
    //! (i.e. in the interval (M_initialTime,M_periodTime)).
    //! Useful for periodic behavior, e.g. BC.
    /*!
     * @return time value in first cycle
     */
    Real inCycleTime() const { return (M_time - static_cast<int>(floor((M_time-M_timeStep/2)/M_periodTime)) * M_periodTime); }

    //! Get the previous time
    /*!
     * @return previous time value
     */
    Real previousTime() const { return M_time - M_timeStep; }

    //! Get the next time
    /*!
     * @return next time value
     */
    Real nextTime() const { return M_time + M_timeStep; }

    //! Get the time step used for advancing
    /*!
     * @return time step value
     */
    const Real& timeStep() const { return M_timeStep; }

    //! Get the number of time step performed
    /*!
     * @return time step performed
     */
    const UInt& timeStepNumber() const { return M_timeStepNumber; }

  //@}

private:

    Real round( const Real& n, const Int& decimal=10 ) const;

    //! initial time
    Real                    M_initialTime;

    //! end time
    Real                    M_endTime;

    //! period time
    Real                    M_periodTime;

    //! current time
    Real                    M_time;

    //! time step
    Real                    M_timeStep;

    //! iteration number
    UInt                    M_timeStepNumber;
};

} // namespace LifeV

#endif // TimeData_H
