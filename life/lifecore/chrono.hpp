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
  @brief Chronometer and fake chronometer class class

  @date 11-12-2010
  @author

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef CHRONO_H
#define CHRONO_H

#include <ctime>

#include <life/lifecore/life.hpp>

namespace LifeV
{

//! @name Chrono - chronometer class
/*!
  This class is used for timing sections of code
*/
class Chrono
{
public:

    //! @name Constructor
    //@{
    Chrono():
        M_t1(0),
        M_t2(0),
        M_dt(0),
        M_running(false)
    {
    }
    //@}

    //! @name Public methods
    //@{
    void reset()
    {
        M_dt = 0;
        M_running = false;
    }

    //! Start the timer
    void start()
    {
        M_t1 = clock();
        M_running = true;
    }

    //! Stop the timer
    void stop()
    {
        if (M_running)
        {
            M_t2 = clock();
            M_dt += M_t2 - M_t1;
            M_running = false;
        }
    }

    //! Compute the difference in time between start and stop
    Real diff()
    {
        if (M_running)
            return ( 1. * ( clock() - M_t1 ) ) / CLOCKS_PER_SEC;

        return ( 1. * ( M_t2 - M_t1 ) ) / CLOCKS_PER_SEC;
    }

    //! Return a cumulative time difference
    Real diffCumul()
    {
        if (M_running)
            return ( 1. * ( M_dt +  clock() - M_t1 ) ) / CLOCKS_PER_SEC;

        return ( 1. * M_dt / CLOCKS_PER_SEC );
    }

    // DEPRECATED
    //! Reset the timer
    void __attribute__((__deprecated__)) Reset()
    {
        M_dt = 0;
        M_running = false;
    }
    Real __attribute__((__deprecated__)) diff_cumul()
    {
        if (M_running)
            return ( 1. * ( M_dt +  clock() - M_t1 ) ) / CLOCKS_PER_SEC;

        return ( 1. * M_dt / CLOCKS_PER_SEC );
    }
    //@}

private:
    //! @name Private members
    //@{
    clock_t M_t1, M_t2, M_dt;
    bool    M_running;
    //@}
};

//! @name ChronoFake - dummy chronometer class
/*!
  When the diff and diff_cumul methods are called, always
  returns -1.0
*/
class ChronoFake
{
public:
    //! @name Constructor
    //@{
    ChronoFake()
    {
    }
    //@}

    //! @name Public methods
    //@{
    void reset()
    {
    }
    void start()
    {
    }
    void stop()
    {
    }
    Real diff()
    {
        return -1.;
    };
    Real diffCumul()
    {
        return -1.;
    };

    // // DEPRECATED
    void __attribute__((__deprecated__)) Reset()
    {
    }
    Real __attribute__((__deprecated__)) diff_cumul()
    {
        return -1.;
    };
    //@}
};

}
#endif // CHRONO_H

