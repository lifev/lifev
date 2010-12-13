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
template <typename T>
class singleton
{
public:
    //! @name Public typedefs
    //@{
    typedef T singleton_type;
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
    static instance_type _S_instance;
    static bool _S_destroyed;
};

// ===================================
// Singleton Implementation
// ===================================

template <class T>
typename singleton<T>::instance_type singleton<T>::_S_instance;

template <class T>
bool singleton<T>::_S_destroyed;

template <class T>
inline T& singleton<T>::instance()
{
    if ( !_S_instance )
    {
        makeInstance();
    }
    return *_S_instance;
}

template <class T>
void singleton<T>::makeInstance()
{
    if ( !_S_instance )
    {
        if ( _S_destroyed )
        {
            lifetime_policy::onDeadReference();
            _S_destroyed = false;
        }
        _S_instance = creation_policy::create();
        lifetime_policy::scheduleDestruction( _S_instance, &destroySingleton );
    }
}

template <class T>
void singleton<T>::destroySingleton()
{
    assert( !_S_destroyed );
    creation_policy::destroy( _S_instance );
    _S_instance = 0;
    _S_destroyed = true;
}

} // Namespace LifeV

#endif // SINGLETON_H
