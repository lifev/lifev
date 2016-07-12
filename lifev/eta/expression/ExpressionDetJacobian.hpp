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
    @brief File containing the Expression for the DETERMINANT JACOBIAN

    @author DAVIDE FORTI
    @date 2016
 */

#ifndef EXPRESSION_DETJACOBIAN_HPP
#define EXPRESSION_DETJACOBIAN_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! ExpressionDetJacobian - Expression for the diameter of the element
/*!
    @author Samuel Quinodoz
 */
class ExpressionDetJacobian : public ExpressionBase< ExpressionDetJacobian >
{
public:

    //! @name Public Types
    //@{

    typedef ExpressionBase<ExpressionDetJacobian> base_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    ExpressionDetJacobian();

    //! Copy constructor
    ExpressionDetJacobian ( const ExpressionDetJacobian& );

    //! Destructor
    ~ExpressionDetJacobian();

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout);

    //@}
};

//! Instance to be used in the expressions
const ExpressionDetJacobian detJ;

} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif /* EXPRESSION_DETJACOBIAN_HPP */
