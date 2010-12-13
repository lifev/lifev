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
  @brief Policies for object lifetime management

  @date 10-09-2004
  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef POLICY_H
#define POLICY_H 1

namespace LifeV
{
/**
   @class policyCreationUsingNew
*/
template <class T>
struct policyCreationUsingNew
{
    //! @name Static methods
    //@{
    static T* create()
    {
        return new T;
    }

    static void destroy( T* p )
    {
        delete p;
    }
    //@}
};

/**
   @class policyLifeTimeDefault
*/
template <class T>
struct policyLifeTimeDefault
{
    //! @name Static methods
    //@{
    static void scheduleDestruction( T*, void ( *pFun ) () )
    {
        std::atexit( pFun );
    }

    static void onDeadReference()
    {
        throw std::logic_error( "Dead Reference Detected" );
    }
    //@}
};

} // Namespace LifeV

#endif // POLICY_H
