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
    @brief The implementation of the h_K expression

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 03 Jan 2012
 */

#include <lifev/eta/expression/ExpressionDetJacobian.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{


// ===================================================
// Constructors & Destructor
// ===================================================
ExpressionDetJacobian::ExpressionDetJacobian()
    : base_Type()
{}

ExpressionDetJacobian::ExpressionDetJacobian ( const ExpressionDetJacobian&)
    : base_Type()
{}

ExpressionDetJacobian::~ExpressionDetJacobian()
{}

// ===================================================
// Methods
// ===================================================

void
ExpressionDetJacobian::display ( std::ostream& out )
{
    out << "detJ";
}


} // Namespace ExpressionAssembly

} // Namespace LifeV
