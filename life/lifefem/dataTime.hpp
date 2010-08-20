//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a class for handling temporal discretization
 *
 *  @author M.A. Fernandez
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 01-06-2009
 */

#ifndef _DATATIME_H_
#define _DATATIME_H_

#include <ostream>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>

namespace LifeV
{


//! DataTime - Class for handling temporal discretization.
/*!
 *  @author Cristiano Malossi
 *
 *  The class is a container for time information.
 */
class DataTime
{
public:

    //! @name Constructors & Destructor
	//@{

    //! Empty Constructor
    DataTime();

    //! Constructor
    /*!
     * @param dataFile GetPot data file
     * @param section the section on the data file that contains the information on the time discretization
     */
    DataTime( const GetPot& dataFile, const std::string& section = "time_discretization" );

    //! Copy constructor
    /*!
     * @param dataTime - DataTime class
     */
    DataTime( const DataTime& dataTime);

    //! Virtual destructor
    virtual ~DataTime() {}

    //@}


    //! @name Methods
    //@{

    //! Update the time by a timestep.
    void updateTime()       { M_time += M_timeStep; }

    //! Return if we can make a new timestep
    /*!
     * @return true if time <= endTime, false viceversa
     */
    bool canAdvance()       { return round( M_time ) <= round( M_endTime ); }

    //! Return if it is the initial time step
    /*!
     * @return true if time = initial time, false viceversa
     */
    bool isFirstTimeStep()  { return round( M_time ) == round( M_initialTime ); }

    //! Return if it is the last time step
    /*!
     * @return true if time + timestep > endTime, false viceversa.
     */
    bool isLastTimeStep()   { return round( M_time + M_timeStep ) > round( M_endTime ); }

    //! Display general information about the content of the class
    /*!
        @param output - specify the output format (std::cout by default)
     */
    virtual void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the initial time step
    /*!
     * @param initialTime initial time step value
     */
    void setInitialTime ( const Real& initialTime ) { M_initialTime = initialTime; }

    //! Set the final time step
    /*!
     * @param endtime final time step value
     */
    void setEndTime     ( const Real& endTime )     { M_endTime = endTime; }

    //! Set the present time of the simulation
    /*!
     * @param time present time value
     */
    void setTime        ( const Real& time )        { M_time = time; }

    //! Set the initial time step
    /*!
     * @param timeStep initial time step value
     */
    void setTimeStep    ( const Real& timeStep )    { M_timeStep = timeStep; }

    //! Set the BDF odert to use
    /*!
     * @param order BDF order
     */
    void setBDF_order   ( const UInt& BDF_order )   { M_BDF_order = BDF_order; }

    //@}


    //! @name Get Methods
	//@{

    //! Get the initial time step
    /*!
     * @return initial time step value
     */
    Real getInitialTime()	 const { return M_initialTime; }

    //! Get the final time step
    /*!
     * @return final time step value
     */
    Real getEndTime()		 const { return M_endTime; }

    //! Get the present time
    /*!
     * @return time value
     */
    Real getTime()			 const { return M_time; }

    //! Get the previous time
    /*!
     * @return previous time value
     */
    Real getPreviousTime()   const { return M_time - M_timeStep; }

    //! Get the next time
    /*!
     * @return next time value
     */
    Real getNextTime()   const { return M_time + M_timeStep; }

    //! Get the time step used for advancing
    /*!
     * @return time step value
     */
    Real getTimeStep()		 const { return M_timeStep; }

    //! Get the number of time step performed
    /*!
     * @return time step performed
     */
    UInt getTimeStepNumber() const { return static_cast <UInt> ( round( round( M_time-M_initialTime ) / round( M_timeStep ) ) ); }

    //! Get the BDF order used
    /*!
     * @return BDF order value
     */
    UInt getBDF_order()		 const { return M_BDF_order; }

    //@}

private:

    Real round( const Real n, const Int decimal=10 ) const;

	Real					M_initialTime;	// initial time
    Real					M_endTime;		// end time
    Real					M_time;			// time
    Real					M_timeStep; 	// time step
    UInt					M_BDF_order; 	// order of the time discretization formula

};

} // namespace LifeV

#endif
