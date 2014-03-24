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
    @brief File containing the Expression for the normal to a surface

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 03 Jan 2012
 */

#ifndef EXPRESSION_NORMAL_HPP
#define EXPRESSION_NORMAL_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! ExpressionNormal - Expression for the normal
/*!
    @author Samuel Quinodoz
 */
class ExpressionNormal : public ExpressionBase< ExpressionNormal >
{
public:

    //! @name Public Types
    //@{

    typedef ExpressionBase<ExpressionNormal> base_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    ExpressionNormal();

    //! Copy constructor
    ExpressionNormal ( const ExpressionNormal& );

    //! Destructor
    ~ExpressionNormal();

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout);

    //@}
};

//! Instance to be used in the expressions
const ExpressionNormal Nface;

} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif /* EXPRESSION_NORMAL_HPP */
