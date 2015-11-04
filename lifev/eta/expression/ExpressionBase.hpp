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
    @brief File where the base class for the expressions is defined

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_BASE_HPP
#define EXPRESSION_BASE_HPP

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionBase  Base class (static polymorphism, CRTP sense) for all the expressions used in assembly procedures.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is nothing else than an empty mother class for all the expressions, in order to
  distinguish those classes from the other ones (acts like a type marker).

  <b> Template parameters </b>

  <i>DerivedType</i>: The type of the derived class.

  <b> Template requirements </b>

  <i>DerivedType</i>: Copiable

*/
template <typename DerivedType>
class ExpressionBase
{
public:

    //! @name Public Types
    //@{

    typedef DerivedType derived_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty and only constructor
    ExpressionBase() {}

    //! Destructor
    virtual ~ExpressionBase() {}

    //@}


    //! @name Methods
    //@{

    //! Method to cast away the type and get the real (DerivedType) object
    const derived_Type& cast() const
    {
        return static_cast<const derived_Type&> (*this);
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No copy (avoid slicing)
    ExpressionBase (const ExpressionBase<DerivedType>&);

    //! No equality (avoid slicing)
    void operator= (const ExpressionBase<DerivedType>&);

    //@}
};


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
