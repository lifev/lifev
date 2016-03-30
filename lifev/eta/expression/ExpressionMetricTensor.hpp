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
    @brief File containing the Expression for the metric tensor of the element

    @author Davide Forti
 */

#ifndef EXPRESSION_METRICTENSOR_HPP
#define EXPRESSION_METRICTENSOR_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

class ExpressionMetricTensor : public ExpressionBase< ExpressionMetricTensor >
{
public:

    //! @name Public Types
    //@{

    typedef ExpressionBase<ExpressionMetricTensor> base_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    ExpressionMetricTensor();

    //! Copy constructor
    ExpressionMetricTensor ( const ExpressionMetricTensor& );

    //! Destructor
    ~ExpressionMetricTensor();

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout);

    //@}
};

//! Instance to be used in the expressions
const ExpressionMetricTensor G;

} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif /* EXPRESSION_METRICTENSOR_HPP */
