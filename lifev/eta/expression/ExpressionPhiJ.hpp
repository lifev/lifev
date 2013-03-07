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
    @brief File containing the definition of the phi_j expression.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_PHI_J_HPP
#define EXPRESSION_PHI_J_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionPhiJ  Class representing the value of the solution in an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
*/
class ExpressionPhiJ : public ExpressionBase<ExpressionPhiJ>
{
public:

    //! @name Public Types
    //@{

    typedef ExpressionBase<ExpressionPhiJ> base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    ExpressionPhiJ();

    //! Copy constructor
    ExpressionPhiJ (const ExpressionPhiJ&);

    //! Destructor
    ~ExpressionPhiJ();

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout);

    //@}
};


//! Simple function to be used in the construction of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
*/
/*inline ExpressionPhiJ
phi_j()
{
    return ExpressionPhiJ();
    };*/

const ExpressionPhiJ phi_j;


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
