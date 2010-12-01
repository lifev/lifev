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
#ifndef _CHRONO_H_INCLUDE
#define _CHRONO_H_INCLUDE

#include <ctime>

namespace LifeV
{
class Chrono
{
public:
    Chrono()
            :
            _t1(0), _t2(0), _dt(0),
            running(false)
    {
    }
    void Reset()
    {
        _dt = 0;
        running = false;
    }
    void start()
    {
        _t1 = clock();
        running = true;
    };
    void stop()
    {
        if (running)
        {
            _t2 = clock();
            _dt += _t2 - _t1;
            running = false;
        }
    };
    double diff()
    {
        if (running)
            return ( 1. * ( clock() - _t1 ) ) / CLOCKS_PER_SEC;

        return ( 1. * ( _t2 - _t1 ) ) / CLOCKS_PER_SEC;
    };
    double diff_cumul()
    {
        if (running)
            return ( 1. * ( _dt +  clock() - _t1 ) ) / CLOCKS_PER_SEC;

        return ( 1. * _dt / CLOCKS_PER_SEC );
    };

private:
    clock_t _t1, _t2, _dt;
    bool    running;

};

class ChronoFake
{
public:
    ChronoFake()
    {
    }
    void Reset()
    {
    }
    void start()
    {
    }
    void stop()
    {
    }
    double diff()
    {
        return -1.;
    };
    double diff_cumul()
    {
        return -1.;
    };

};

}
#endif

