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
  @brief Wall clock timer class

  @date 01.10.2010
  @author Radu Popescu <radu.popescu@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef WALL_CLOCK_H
#define WALL_CLOCK_H 1

#include <sys/time.h>

namespace LifeV
{

//! Wall clock timer class
/*!
  ... totul e minunat ...
*/
class WallClock
{
public:
    //! @name Public typedefs
    //@{
    typedef struct timeval time_Type;
    //@}

    //! @name Constructor and destructor
    //@{
    WallClock();
    ~WallClock();
    //@}

    //! @name Public methods
    //@{
    void start();
    void stop();
    void reset();
    //@}

    //! @name Get methods
    //@{
    const double& elapsedTime() const
    {
        return M_elapsedTime;
    }
    //@}

private:
    // Disabled copy constructor and assignment operator
    WallClock (const WallClock&);
    WallClock& operator= (const WallClock&);

    time_Type M_startTime;
    time_Type M_stopTime;

    double M_elapsedTime;
};

} // namespacec LifeV

#endif // WALL_CLOCK_H
