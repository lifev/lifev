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
    @brief

    @author Antonio Cervone <ant.cervone@gmail.com>
    @date 2013-07-15
 */
#ifndef LIFECHRONOMANAGER_HPP
#define LIFECHRONOMANAGER_HPP

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChrono.hpp>

namespace LifeV
{

template <typename TimerType = LifeChrono>
class LifeChronoManager
{
public:
    typedef TimerType timer_Type;
    typedef std::map<std::string const, timer_Type*> timerList_Type;
    typedef boost::shared_ptr<Epetra_Comm const> commPtr_Type;

    LifeChronoManager( commPtr_Type comm ):
        M_stringMaxSize( 16 ),
        M_comm( comm )
    {}

    void add( std::string const & name, timer_Type* timer )
    {
        UInt const nameSize = name.size();
        if( nameSize > M_stringMaxSize ) M_stringMaxSize = nameSize;
        M_timerList.insert( std::make_pair( name, timer) );
    }

    void print( std::ostream & out = std::cout )
    {
        UInt const printSize = 80;
        bool isLeader = M_comm->MyPID() == 0;

        std::vector<Real> times( M_timerList.size() );
        Real globalTime = 0;
        UInt count = 0;
        for( typename timerList_Type::const_iterator it = M_timerList.begin();
             it != M_timerList.end(); ++it, count++ )
        {
            times[ count ] = it->second->globalDiff( *M_comm );
            globalTime += times[ count ];
        }

        if( isLeader )
        {
            out << std::string(printSize, '=') << std::endl;
            out << std::setw( M_stringMaxSize ) << "Name";
            out << std::setw( 16 ) << "Time (s)";
            out << std::setw( 16 ) << "Perc (%)" << std::endl;
            out << std::string(printSize, '=') << std::endl;

            count = 0;
            for( typename timerList_Type::const_iterator it = M_timerList.begin();
                 it != M_timerList.end(); ++it, count++ )
            {
                out << std::setw( M_stringMaxSize ) << it->first;
                out << std::setw( 16 ) << std::fixed << std::setprecision(2) << times[ count ];
                out << std::setw( 16 ) << std::fixed << std::setprecision(2) << 100.* times[ count ] / globalTime << std::endl;
            }
            out << std::string(printSize, '=') << std::endl;
            out << std::setw( M_stringMaxSize ) << "Total Time";
            out << std::setw( 16 ) << std::fixed << std::setprecision(2) << globalTime;
            out << std::setw( 16 ) << std::fixed << std::setprecision(2) << 100. << std::endl;
            out << std::string(printSize, '=') << std::endl;
        }
    }

protected:
    timerList_Type M_timerList;
    UInt M_stringMaxSize;
    commPtr_Type M_comm;
};

} // end namespace LifeV

#endif // LIFECHRONOMANAGER_HPP
