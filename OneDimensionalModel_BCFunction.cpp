//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing the interface class for the boundary function of 1D model.
 *
 *  @version 1.0
 *  @author Lucia Mirabella
 *  @date 01-08-2006
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-04-2010
 *
 */

#include <lifemc/lifefem/OneDimensionalModel_BCFunction.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCFunction::OneDimensionalModel_BCFunction() :
    M_function  ()
{}

OneDimensionalModel_BCFunction::OneDimensionalModel_BCFunction( const Function_Type& function ) :
    M_function  ( function )
{}

OneDimensionalModel_BCFunction::OneDimensionalModel_BCFunction( const OneDimensionalModel_BCFunction& BCFunction ) :
    M_function  ( BCFunction.M_function )
{}

// ===================================================
// Operators
// ===================================================
OneDimensionalModel_BCFunction&
OneDimensionalModel_BCFunction::operator=( const OneDimensionalModel_BCFunction& BCFunction )
{
    if ( this != &BCFunction )
    {
        M_function = BCFunction.M_function;
    }

    return *this;
}

Real
OneDimensionalModel_BCFunction::operator()( const Real& time, const Real& timeStep ) const
{
    return M_function( time, timeStep );
}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_BCFunction::setFunction( const Function_Type& function )
{
    M_function = function;
}

// ===================================================
// Get Methods
// ===================================================
const OneDimensionalModel_BCFunction::Function_Type&
OneDimensionalModel_BCFunction::Function() const
{
    return M_function;
}

}
