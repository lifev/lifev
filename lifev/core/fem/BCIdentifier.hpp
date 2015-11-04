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

#ifndef BCIDENTIFIER_H
#define BCIDENTIFIER_H 1

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

//! BCIdentifierBase - Base class holding DOF identifiers for implementing BC

class BCIdentifierBase
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    BCIdentifierBase()
    {
        // Nothing to be done here
    }

    //! Constructor given the ID
    /*!
        Creates an BCIdentifier with a given ID
        @param i Usually the id of the DOF, or the id of a boundary face, etc...
    */
    explicit BCIdentifierBase ( ID const& i ) : M_id ( i )
    {
        // Nothing to be done here
    }

    //! Copy constructor
    BCIdentifierBase ( BCIdentifierBase const& id );

    //! Destructor
    virtual ~BCIdentifierBase()
    {
        // Nothing to be done here
    }

    //@}

    //! @name Methods
    //@{

    //! Display method
    virtual void showMe ( std::ostream& output = std::cout ) const;

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

    //! Returns the ID of the BCIdentifier
    const ID& id() const
    {
        return M_id;
    }

    //@}

protected:

    ID M_id;
};

//! BCIdentifierComp - Functor for ordering operations (required in set STL container)

class BCIdentifierComparison
{
public:

    //! @name Operators
    //@{

    //! Comparison operator for BCIdentifier objects
    /*!
        @return Boolean which is true if the ID of the first BCIdentifier is smaller
                than the ID of the second BCIdentifier
     */
    bool operator() ( const BCIdentifierBase* i1, const BCIdentifierBase* i2 ) const
    {
        return ( i1->id() < i2->id() );
    }

    //! Comparison operator for shared pointers to BCIdentifier objects
    /*!
        @return Boolean which is true if the ID of the first BCIdentifier is smaller
                than the ID of the second BCIdentifier
     */
    bool operator() ( std::shared_ptr<BCIdentifierBase> const& i1, std::shared_ptr<BCIdentifierBase> const& i2 ) const
    {
        return ( i1.get()->id() < i2.get()->id() );
    }

    //@}
};

//! Overloading == operator for objects of type BCIdentifier
/*!
    @param first The first BCIdentifier
    @param second The second BCIdentifier
    @return A bool which is 1 if the ID of the two BCIdentifier objects are the same
 */
inline bool operator== ( const BCIdentifierBase& first, const BCIdentifierBase& second )
{
    return first.id() == second.id();
}

//! BCIdentifierEssential - BCIdentifier for implementing Essential Boundary Conditions
/*!

    This class holds the DOF identifier and its coordinates for implementing Essential
    Boundary Conditions

 */

class BCIdentifierEssential: public BCIdentifierBase
{
public:

    //! Constructor & Destructor
    //@{

    //! Empty Constructor
    BCIdentifierEssential() : BCIdentifierBase()
    {
        // Nothing to be done here
    }

    //! Constructor given the ID and the coordinates
    /*!
     *  Creates an BCIdentifier with a given ID and given coordinates
        @paramx x x-coordinate of the node where this BC applies
        @paramx y y-coordinate of the node where this BC applies
        @paramx z z-coordinate of the node where this BC applies
     */
    BCIdentifierEssential ( const ID& id, const Real& x, const Real& y, const Real& z ) :
        BCIdentifierBase ( id ),
        M_x ( x ),
        M_y ( y ),
        M_z ( z )
    {
        // Nothing to be done here
    }

    //! Copy Constructor
    BCIdentifierEssential ( BCIdentifierEssential const& id );

    //@}

    //! @name Methods
    //@{

    //! Display method
    virtual void showMe ( std::ostream& output = std::cout ) const;

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


//! BCIdentifierNatural - Idenifier for Natural and Robin Boundary Condiions
/*!

    This class holds the DOF identifier and the bdLocalToGlobal information for implementing
    Natural and Robin boundary conditions

 */

class BCIdentifierNatural: public BCIdentifierBase
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    BCIdentifierNatural() : BCIdentifierBase()
    {
        // Nothing to be done here
    }

    //! Constructor given ID and bdLocalToGlobal map
    /*!
        Creates an BCIdentifier with a given ID and a given local-to-global map
        @param i The number of the boundary face
        @param localToGlobal A vector holding the local-to-global map on this face
    */
    BCIdentifierNatural ( const ID& i, const std::vector<ID>& localToGlobal );

    //! Constructor given the ID
    /*!
        @param id The ID of the dof
    */
    explicit BCIdentifierNatural ( const ID& id );

    //! Copy Constructor
    BCIdentifierNatural ( BCIdentifierNatural const& id );

    //! Destructor
    virtual ~BCIdentifierNatural()
    {
        // Nothing to be done here
    }

    //@}

    //! @name Methods
    //@{

    //! Display method
    virtual void showMe (std::ostream& output = std::cout ) const;

    //@}

    //! @name Get Methods
    //@{

    //! Return the global DOF corresponding tho the i-th local DOF in the face
    /*!
        @param i The local DOF in the face
    */
    ID boundaryLocalToGlobalMap ( const ID& i ) const
    {
        return M_localToGlobal [i ];
    }

    //@}

private:

    std::vector<ID> M_localToGlobal;
};

} // Namespace LifeV

#endif /* BCIDENTIFIER_H */
