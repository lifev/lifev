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
    @brief File containing the definition of the lapacian(phi_i) expression.

    @author Davide Forti <davide.forti@epfl.ch>

    @date 11-2014
 */

#ifndef EXPRESSION_LAPLACIAN_J_HPP
#define EXPRESSION_LAPLACIAN_J_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>
#include <lifev/eta/expression/ExpressionPhiJ.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionLaplacianJ  Class representing the laplacian of the trial function in an expression
/*!
  @author Davide Forti <davide.forti@epfl.ch>
*/
class ExpressionLaplacianJ : public ExpressionBase<ExpressionLaplacianJ>
{
public:

    //! @name Public Types
    //@{

    typedef ExpressionBase<ExpressionLaplacianJ> base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    ExpressionLaplacianJ();

    //! Copy constructor
    ExpressionLaplacianJ (const ExpressionLaplacianJ&);

    //! Destructor
    ~ExpressionLaplacianJ();

    //@}

    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout);

    //@}
};

//! Simple function to be used in the construction of an expression
/*!
  @author Davide Forti <davide.forti@epfl.ch>
*/
inline ExpressionLaplacianJ
laplacian (const ExpressionPhiJ& )
{
    return ExpressionLaplacianJ();
}


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
