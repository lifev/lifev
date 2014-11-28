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
    @brief Class representing the laplacian of the test function in an expression

    @author Davide Forti <davide.forti@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_LAPLACIAN_I_CPP
#define EXPRESSION_LAPLACIAN_I_CPP


#include <lifev/eta/expression/ExpressionLaplacianI.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

// ===================================================
// Constructors & Destructor
// ===================================================

ExpressionLaplacianI::ExpressionLaplacianI()
    : base_Type()
{}


ExpressionLaplacianI::ExpressionLaplacianI (const ExpressionLaplacianI&)
    : base_Type()
{}


ExpressionLaplacianI::~ExpressionLaplacianI()
{}

// ===================================================
// Methods
// ===================================================

void
ExpressionLaplacianI::display (std::ostream& out)
{
    out << "dphi_i";
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
