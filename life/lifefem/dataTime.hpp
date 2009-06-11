/*
This file is part of the LifeV library
Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*!
  \file dataTime.hpp

  \version 1.0
  \date 01/2003
  \author M.A. Fernandez

  \version 1.9
  \date 06/2009
  \author Cristiano Malossi<cristiano.malossi@epfl.ch>

  \brief File containing a class for handling temporal discretization
*/

#ifndef _DATATIME_H_
#define _DATATIME_H_

#include <ostream>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>

namespace LifeV
{

/*!
 * \class DataTime
 * \brief Base class which holds data concerning temporal discretization
 *
 * @author M.A. Fernandez, Cristiano Malossi
 * @see
 */
class DataTime
{
public:

    /** @name Constructors, Destructor
     */
	//@{

    DataTime( const GetPot& dfile, const std::string& section = "time" );

    DataTime( const DataTime& dataTime);

    //! Virtual destructor
    virtual ~DataTime() {}

    //@}



    /** @name Get Functions
     */
	//@{

    Real getInitialTime()	const { return M_initialTime; }

    Real getEndTime()		const { return M_endTime; }

    Real getTime()			const { return M_time; }

    Real getTimeStep()		const { return M_timeStep; }

    UInt getBDF_order()		const { return M_BDF_order; }

    //@}



    /** @name Set Functions
     */
	//@{

    void setInitialTime	( const Real& initialTime ) { M_initialTime = initialTime; }

    void setEndTime		( const Real& endTime ) 	{ M_endTime = endTime; }

    void setTime		( const Real& time ) 		{ M_time = time; }

    void setTimeStep	( const Real& timeStep ) 	{ M_timeStep = timeStep; }

    void setBDF_order	( const UInt& BDF_order ) 	{ M_BDF_order = BDF_order; }

    //@}



    //! Ouptut
    virtual void showMe( std::ostream& output = std::cout ) const;

private:

	Real					M_initialTime;	// initial time
    Real					M_endTime;		// end time
    Real					M_time;			// time
    Real					M_timeStep; 	// time step
    UInt					M_BDF_order; 	// order of the time discretization formula

};

} // namespace LifeV

#endif
