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

#include "WallClock.hpp"

#include <cstdlib>

namespace LifeV
{
// ==========================
// Constructor and destructor
// ==========================
WallClock::WallClock()
{
    reset();
}

WallClock::~WallClock()
{
}

// ==========================
// Public methods
// ==========================

void WallClock::start()
{
    gettimeofday (&M_startTime, NULL);
}

void WallClock::stop()
{
    gettimeofday (&M_stopTime, NULL);

    time_t seconds = M_stopTime.tv_sec - M_startTime.tv_sec;
    suseconds_t microseconds = M_stopTime.tv_usec - M_startTime.tv_usec;

    M_elapsedTime += static_cast<double> (seconds + microseconds / 1000000.0);
}

void WallClock::reset()
{
    M_elapsedTime = 0;
}

} // namespace LifeV
