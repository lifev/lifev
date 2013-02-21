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

#ifndef FACTORY_SINGLETON_H
#define FACTORY_SINGLETON_H 1

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <exception>
#include <new>
#include <stdexcept>

#include <lifev/core/util/FactoryPolicy.hpp>

namespace LifeV
{
/**
   \class FactorySingleton
   \brief implement the FactorySingleton pattern

   A FactorySingleton pattern implementation using the ideas
   from Alexandrescu's book "modern C++ design"
   http://www.moderncppdesign.com/
*/
template <typename SingletonType>
class FactorySingleton
{
public:
    //! @name Public typedefs
    //@{
    typedef SingletonType Singleton_Type;
    typedef FactoryPolicyLifeTimeDefault<Singleton_Type> lifetimeFactoryPolicy_Type;
    typedef FactoryPolicyCreationUsingNew<Singleton_Type> creationFactoryPolicy_Type;
    //@}

    //! @name Public static methods
    //@{
    /**
       return the instance of the FactorySingleton
    */
    static Singleton_Type& instance();
    //@}
private:
    //! @name Private typedefs
    //@{
    typedef Singleton_Type* instance_Type;
    //@}

    // Disable the constructor
    FactorySingleton();

    //! @name Private static methods
    //@{
    static void makeInstance();

    /**
       FactorySingleton::makeInstance (helper for Instance)
    */
    static void destroyFactorySingleton();

    //@}
    static instance_Type S_instance;
    static bool S_destroyed;
};

// ===================================
// FactorySingleton Implementation
// ===================================

template <class SingletonType>
typename FactorySingleton<SingletonType>::instance_Type FactorySingleton<SingletonType>::S_instance;

template <class SingletonType>
bool FactorySingleton<SingletonType>::S_destroyed;

template <class SingletonType>
inline SingletonType& FactorySingleton<SingletonType>::instance()
{
    if ( !S_instance )
    {
        makeInstance();
    }
    return *S_instance;
}

template <class SingletonType>
void FactorySingleton<SingletonType>::makeInstance()
{
    if ( !S_instance )
    {
        if ( S_destroyed )
        {
            lifetimeFactoryPolicy_Type::onDeadReference();
            S_destroyed = false;
        }
        S_instance = creationFactoryPolicy_Type::create();
        lifetimeFactoryPolicy_Type::scheduleDestruction ( S_instance, &destroyFactorySingleton );
    }
}

template <class SingletonType>
void FactorySingleton<SingletonType>::destroyFactorySingleton()
{
    assert ( !S_destroyed );
    creationFactoryPolicy_Type::destroy ( S_instance );
    S_instance = 0;
    S_destroyed = true;
}

} // Namespace LifeV

#endif // FACTORY_SINGLETON_H
