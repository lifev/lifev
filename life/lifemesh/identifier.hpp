/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
  \file identifier.h
  \brief Classes for identifiers
  \version 1.0
  \author M.A. Fernandez & Luca Formaggia
  \date 07/2002

  This classes hold a identifier that allow us to impose a specific boundary condition.
  Each type of boundary condition needs a specic information on the boundary. Thus, the
  key is to use inheritance by adding, to the base class, the information requested for
  imposing the BC.

*/

#ifndef __IDENTIFIER_HH__
#define __IDENTIFIER_HH__

#include <iostream>

#include <boost/shared_ptr.hpp>

#include "lifeV.hpp"
#include "SimpleVect.hpp"

namespace LifeV
{

//============ IdentifierBase ==============

/*!

 \class IdentifierBase

 Base class holding Dof identifiers for implementing BC

 \todo The data functions given by the user must have the following declaration
 \verbatim
 Real g(const Real& time, const Real& x, const Real& y, const Real& z, const ID& icomp)
 \endverbatim

 We can use inheritance to hold specific boundary condition data. See, for instance,
 Mixte boundary conditions.
*/
//! Declaration of the base class holding DOF identifiers for implementing BC
class IdentifierBase
{
public:

    //! Constructor
    /*!
      \param i ussualy the id of the Dof, id of a boundary face, etc...
    */
    IdentifierBase( ID const & i ) : _id( i )
        {
            //nothing to be done here
        }

    virtual ~IdentifierBase()
        {
            //nothing to be done here
        }

    //! Returns the ID
    const ID& id() const
    {
        return _id;
    }

    //! Conversion operators
    operator unsigned int()
    {
        return id();
    }
    operator int()
    {
        return ( int ) id();
    }

protected:
    //! The identifier
    ID _id;
};

/*!

 \class identifierComp

 Functor for ordering operations (requested in set STL container)

*/
class identifierComp
{
public:
    bool operator() ( const IdentifierBase* i1, const IdentifierBase* i2 ) const
    {
        return ( i1->id() < i2->id() );
    }
    bool operator() ( boost::shared_ptr<IdentifierBase> const& i1, boost::shared_ptr<IdentifierBase> const& i2 ) const
    {
        return ( i1.get()->id() < i2.get()->id() );
    }
};

//! Overloading == operator for identifiers
inline bool operator==( const IdentifierBase& a, const IdentifierBase& b )
{
    return a.id() == b.id();
}


//============ IdentifierEssential ==============
/*!

 \class IdentifierEssential

 Class holding the Dof identifier and coordinates for implementing Essential BC

 \todo A Essential boundary condition requests the number of the Dof (in a scalar sense) and the
 coordiantes of associated node. This information is updated in Dof::bdUpdate method
*/
class IdentifierEssential: public IdentifierBase
{
public:

    //! Constructor
    /*!
      \param i the id of the Dof
      \param x x-coordinate of the node where this BC applies
      \param y y-coordinate of the node where this BC applies
      \param z z-coordinate of the node where this BC applies
    */
    IdentifierEssential( const ID& id, const Real& x, const Real& y, const Real& z ) : IdentifierBase( id )
    {
        _x = x;
        _y = y;
        _z = z;
    }

    //! Recovering node coordinates
    const Real& x() const
    {
        return _x;
    }
    const Real& y() const
    {
        return _y;
    }
    const Real& z() const
    {
        return _z;
    }

private:
    //! Node coordinates
    Real _x, _y, _z;
};


//============ IdentifierNatural ==============

/*!

 \class IdentifierNatural

 Class holding the Dof identifier and the bdLocalToGlobal information for implementing
 Natural and Mixte boundary conditions

 \todo Natural or Mixte boundary conditions requests the number of the boundary face number
 where they apply and the bdLocalToGlobal map on this face. This information is updated
 in Dof::bdUpdate method
*/
class IdentifierNatural: public IdentifierBase
{
public:

    //! Constructor
    /*!
      \param i the number of the boundary face
      \param bdltg a SimpleVect holding the bdLocalToGlobal map on this face
    */
    IdentifierNatural( const ID& i, const SimpleVect<ID>& bdltg );

    //! Constructor when a vector data is provided
    /*!
      \param i the number of the dof
    */
    IdentifierNatural( const ID& i );


    //! Return the global Dof corresponding tho the i-th local Dof in the face
    /*!
      \param i local Dof in the face
    */
    ID bdLocalToGlobal( const ID& i ) const
    {
        return _bdltg( i );
    }

private:
    //! SimpleVect container holding the bdLocalToGlobal map on this face
    SimpleVect<ID> _bdltg;
};

}
#endif
