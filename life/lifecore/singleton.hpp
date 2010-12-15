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
  @brief Singleton pattern class

  @date 10-09-2004
  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef SINGLETON_H
#define SINGLETON_H 1

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <exception>
#include <new>
#include <stdexcept>

#include <life/lifecore/policy.hpp>

namespace LifeV
{
/**
   \class singleton
   \brief implement the Singleton pattern

   A Singleton pattern implementation using the ideas
   from Alexandrescu's book "modern C++ design"
   http://www.moderncppdesign.com/
*/
template <typename SingletonType>
class singleton
{
public:
    //! @name Public typedefs
    //@{
    typedef SingletonType singleton_type;
    typedef policyLifeTimeDefault<singleton_type> lifetime_policy;
    typedef policyCreationUsingNew<singleton_type> creation_policy;
    //@}

    //! @name Public static methods
    //@{
    /**
       return the instance of the singleton
    */
    static singleton_type& instance();
    //@}
private:
    //! @name Private typedefs
    //@{
    typedef singleton_type* instance_type;
    //@}

    // Disable the constructor
    singleton();

    //! @name Private static methods
    //@{
    static void makeInstance();

    /**
       Singleton::makeInstance (helper for Instance)
    */
    static void destroySingleton();

    //@}
    static instance_type S_instance;
    static bool S_destroyed;
};

// ===================================
// Singleton Implementation
// ===================================

template <class SingletonType>
typename singleton<SingletonType>::instance_type singleton<SingletonType>::S_instance;

template <class SingletonType>
bool singleton<SingletonType>::S_destroyed;

template <class SingletonType>
inline SingletonType& singleton<SingletonType>::instance()
{
    if ( !S_instance )
    {
        makeInstance();
    }
    return *S_instance;
}

template <class SingletonType>
void singleton<SingletonType>::makeInstance()
{
    if ( !S_instance )
    {
        if ( S_destroyed )
        {
            lifetime_policy::onDeadReference();
            S_destroyed = false;
        }
        S_instance = creation_policy::create();
        lifetime_policy::scheduleDestruction( S_instance, &destroySingleton );
    }
}

template <class SingletonType>
void singleton<SingletonType>::destroySingleton()
{
    assert( !S_destroyed );
    creation_policy::destroy( S_instance );
    S_instance = 0;
    S_destroyed = true;
}

} // Namespace LifeV

#endif // SINGLETON_H
