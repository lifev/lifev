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
    @brief File containing the Expression for the position

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 03 Jan 2012
 */

#ifndef EXPRESSION_POSITION_HPP
#define EXPRESSION_POSITION_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! ExpressionPosition - Expression for the position (=quadrature node)
/*!
    @author Samuel Quinodoz
 */
class ExpressionPosition : public ExpressionBase< ExpressionPosition >
{
public:

    //! @name Public Types
    //@{

    typedef ExpressionBase<ExpressionPosition> base_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    ExpressionPosition();

    //! Copy constructor
    ExpressionPosition ( const ExpressionPosition& );

    //! Destructor
    ~ExpressionPosition();

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout);

    //@}
};

//! Instance to be used in the expressions
const ExpressionPosition X;

} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif /* EXPRESSION_POSITION_HPP */
