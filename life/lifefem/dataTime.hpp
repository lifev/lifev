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

#ifndef _DATATIME_H_
#define _DATATIME_H_ 1

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

    //! Return if it is the last time step
    /*!
     * @return true if time + timestep > endTime, false viceversa.
     */
    bool isLastTimeStep() { return round( M_time + M_timeStep ) > round( M_endTime ); }

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

    //! Set the BDF odert to use
    /*!
     * @param order BDF order
     */
    void setOrderBDF( const UInt& orderBDF ) { M_orderBDF = orderBDF; }
    void __attribute__ ((__deprecated__)) setBDF_order( const UInt& orderBDF ) { setOrderBDF(orderBDF); }

    //! Set the theta of Newmark scheme
    /*!
     * @param theta - coefficient of Newmark scheme
     */
    void setTheta( const Real& theta ) { M_theta = theta; }

    //! Set the theta of Newmark scheme
    /*!
     * @param gamma- coefficient of Newmark scheme
     */
    void setGamma( const Real& gamma ) { M_gamma = gamma; }
    void __attribute__ ((__deprecated__)) setZeta ( const Real& gamma ) { setGamma(gamma); }

    //@}


    //! @name Get Methods
    //@{

    //! Get the initial time step
    /*!
     * @return initial time step value
     */
    const Real& initialTime() const { return M_initialTime; }
    const Real& __attribute__ ((__deprecated__)) getInitialTime() const { return initialTime(); }

    //! Get the final time step
    /*!
     * @return final time step value
     */
    const Real& endTime() const { return M_endTime; }
    const Real& __attribute__ ((__deprecated__)) getEndTime() const{return endTime();}

    //! Get the period
    /*!
     * @return period value
     */
    const Real& periodTime() const { return M_periodTime; }
    const Real& __attribute__ ((__deprecated__))  getPeriodTime() const {return periodTime();}
    //! Get the present time
    /*!
     * @return time value
     */
    const Real& time() const { return M_time; }
    const Real& __attribute__ ((__deprecated__)) getTime() const { return time(); }

    //! Get the time left
    /*!
     * @return time left value
     */
    Real leftTime() const { return round( M_endTime - M_time ); }
    Real getLeftTime() const { return leftTime();}
    //! Get the elapsed time
    /*!
     * @return elapsed time value
    */
    Real elapsedTime() const { return round( M_time - M_initialTime ); }
  Real __attribute__((__deprecated__))  getElapsedTime() const{ return  elapsedTime();}

    //! Get the present time shifted inside the first cycle
    //! (i.e. in the interval (M_initialTime,M_periodTime)).
    //! Useful for periodic behavior, e.g. BC.
    /*!
     * @return time value in first cycle
     */
    Real inCycleTime() const { return (M_time - static_cast<int>(floor((M_time-M_timeStep/2)/M_periodTime)) * M_periodTime); }
    Real __attribute__ ((__deprecated__))getInCycleTime() const { return inCycleTime();}

                    //! Get the previous time
    /*!
     * @return previous time value
     */
    Real previousTime()   const { return M_time - M_timeStep; }
    Real __attribute__ ((__deprecated__)) getPreviousTime()   const { return previousTime(); }

    //! Get the next time
    /*!
     * @return next time value
     */
    Real nextTime()       const { return M_time + M_timeStep; }
    Real __attribute__ ((__deprecated__)) getNextTime() const { return nextTime(); }

    //! Get the time step used for advancing
    /*!
     * @return time step value
     */
    const Real& timeStep()       const { return M_timeStep; }
    const Real& __attribute__ ((__deprecated__))   getTimeStep()		 const { return timeStep(); }

    //! Get the number of time step performed
    /*!
     * @return time step performed
     */
    const UInt& timeStepNumber() const { return M_timeStepNumber; }
    const UInt& __attribute__ ((__deprecated__)) getTimeStepNumber() const { timeStepNumber();}

    //! Get the BDF order used
    /*!
     * @return BDF order value
     */
    const UInt& orderBDF()		 const { return M_orderBDF; }
    const UInt& __attribute__ ((__deprecated__)) getBDF_order()      const { return orderBDF(); }

    //! Return theta parameter of Newmark scheme
    /*!
     * @return theta value
     */
    const Real& theta()             const { return M_theta; }

    //! Return gamma of Newmark scheme
    /*!
     * @return gamma value
     */
    const Real& gamma()              const { return M_gamma; }
  const Real& __attribute__ ((__deprecated__)) zeta()     const { return gamma(); }
    //! Return Newmark parameters (\f$theta\f$, $\gamma$)

    std::vector<Real> coefficientsNewmark();
  std::vector<Real>  __attribute__ ((__deprecated__)) getNewmark_parameters();
    //@}

private:

    Real round( const Real n, const Int decimal=10 ) const;

    Real					M_initialTime;	  // initial time
    Real					M_endTime;		  // end time
    Real					M_periodTime;	  // period time
    Real					M_inCycleTime;	  // in cycle time
    Real					M_time;           // time
    Real					M_timeStep; 	  // time step
    UInt                                     M_timeStepNumber; // iteration number
    UInt					M_orderBDF; 	  // order of the time discretization formula
    Real                                    M_theta;          // Newmark parameter
    Real                                    M_gamma;           // Newmark parameter
};

} // namespace LifeV

#endif
