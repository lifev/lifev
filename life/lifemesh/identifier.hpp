//@HEADER
/*
************************************************************************

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

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Classes for identifiers

    @date 01-07-2002
    @author M. A. Fernandez & Luca Formaggia
    @contributor Luca Bertagna <lbertag@emory.edu>

    This classes hold a identifier that allow us to impose a specific boundary condition.
    Each type of boundary condition needs a specic information on the boundary. Thus, the
    key is to use inheritance by adding, to the base class, the information requested for
    imposing the BC.
 */

#ifndef IDENTIFIER_H
#define IDENTIFIER_H 1

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifearray/SimpleVect.hpp>

namespace LifeV {

//! IdentifierBase - Base class holding Dof identifiers for implementing BC

class IdentifierBase
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    IdentifierBase()
    {
        // Nothing to be done here
    }

    //! Constructor given the ID
    /*!
        Creates an Identifier with a given ID
        @param i Usually the id of the Dof, or the id of a boundary face, etc...
    */
    explicit IdentifierBase( ID const & i ) : M_id( i )
    {
        // Nothing to be done here
    }

    //! Copy constructor
    IdentifierBase( IdentifierBase const & id );

    //! Destructor
    virtual ~IdentifierBase()
    {
        // Nothing to be done here
    }

    //@}

    //! @name Methods
    //@{

    //! Display method
    virtual void showMe( std::ostream & output = std::cout ) const;

    //@}

    //! @name Operators
    //@{

    //! Conversion operator to unsigned int
    operator unsigned int()
    {
        return id();
    }

    //! Conversion operator to int
    operator int()
    {
        return ( int ) id();
    }

    //@}


    //! @name Get Methods
    //@{

    //! Returns the ID of the Identifier
    const ID& id() const
    {
        return M_id;
    }

    //@}

protected:

    ID M_id;
};

//! IdentifierComp - Functor for ordering operations (required in set STL container)

class identifierComp
{
public:

    //! @name Operators
    //@{

    //! Comparison operator for Identifier objects
    /*!
        @return Boolean which is true if the ID of the first Identifier is smaller
                than the ID of the second Identifier
     */
    bool operator() ( const IdentifierBase* i1, const IdentifierBase* i2 ) const
    {
        return ( i1->id() < i2->id() );
    }

    //! Comparison operator for shared pointers to Identifier objects
    /*!
        @return Boolean which is true if the ID of the first Identifier is smaller
                than the ID of the second Identifier
     */
    bool operator() ( boost::shared_ptr<IdentifierBase> const & i1, boost::shared_ptr<IdentifierBase> const & i2 ) const
    {
        return ( i1.get()->id() < i2.get()->id() );
    }

    //@}
};

//! Overloading == operator for objects of type Identifier
/*!
    @param first The first Identifier
    @param second The second Identifier
    @return A bool which is 1 if the ID of the two Identifier objects are the same
 */
inline bool operator==( const IdentifierBase& first, const IdentifierBase& second )
{
    return first.id() == second.id();
}

//! IdentifierEssential - Identifier for implementing Essential Boundary Conditions
/*!

    This class holds the Dof identifier and its coordinates for implementing Essential
    Boundary Conditions

 */

class IdentifierEssential: public IdentifierBase
{
public:

    //! Constructor & Destructor
    //@{

    //! Empty Constructor
    IdentifierEssential() : IdentifierBase()
    {
        // Nothing to be done here
    }

    //! Constructor given the ID and the coordinates
    /*!
     *  Creates an Identifier with a given ID and given coordinates
        @paramx x x-coordinate of the node where this BC applies
        @paramx y y-coordinate of the node where this BC applies
        @paramx z z-coordinate of the node where this BC applies
     */
    IdentifierEssential( const ID& id, const Real& x, const Real& y, const Real& z ) :
            IdentifierBase( id ),
            M_x( x ),
            M_y( y ),
            M_z( z )
    {
        // Nothing to be done here
    }

    //! Copy Constructor
    IdentifierEssential( IdentifierEssential const & id );

    //@}

    //! @name Methods
    //@{

    //! Display method
    virtual void showMe( std::ostream& output = std::cout ) const;

    //@}

    //! @name Get Methods
    //@{

    //! Recovering the node's x-coordinate
    /*!
        @return The x-coordinate of the node
     */
    const Real& x() const
    {
        return M_x;
    }

    //! Recovering the node's y-coordinate
    /*!
        @return The y-coordinate of the node
    */
    const Real& y() const
    {
        return M_y;
    }

    //! Recovering the node's z-coordinate
    /*!
        @return The z-coordinate of the node
     */
    const Real& z() const
    {
        return M_z;
    }

    //@}

private:

    Real M_x, M_y, M_z;
};


//! IdentifierNatural - Idenifier for Natural and Mixte Boundary Condiions
/*!

    This class holds the Dof identifier and the bdLocalToGlobal information for implementing
    Natural and Mixte boundary conditions

 */

class IdentifierNatural: public IdentifierBase
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    IdentifierNatural() : IdentifierBase()
    {
        // Nothing to be done here
    }

    //! Constructor given ID and bdLocalToGlobal map
    /*!
        Creates an Identifier with a given ID and a given local-to-global map
        @param i The number of the boundary face
        @param bdltg A SimpleVect holding the local-to-global map on this face
    */
    IdentifierNatural( const ID& i, const SimpleVect<ID>& localToGlobal );

    //! Constructor given the ID
    /*!
        @param id The ID of the dof
    */
    explicit IdentifierNatural( const ID& id );

    //! Copy Constructor
    IdentifierNatural( IdentifierNatural const & id );

    //! Destructor
    virtual ~IdentifierNatural()
    {
        // Nothing to be done here
    }

    //@}

    //! @name Methods
    //@{

    //! Display method
    virtual void showMe(std::ostream& output = std::cout ) const;

    //@}

    //! @name Get Methods
    //@{

    //! Return the global Dof corresponding tho the i-th local Dof in the face
    /*!
        @param i The local Dof in the face
    */
    ID localToGlobalMap( const ID& i ) const
    {
        return M_localToGlobal( i );
    }

    ID __attribute__ ((__deprecated__)) bdLocalToGlobal( const ID& i ) const
    {
        return localToGlobalMap( i );
    }

    //@}

private:

    SimpleVect<ID> M_localToGlobal;
};

} // Namespace LifeV

#endif /* IDENTIFIER_H */
