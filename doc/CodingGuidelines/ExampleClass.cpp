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
    @brief A short description of the file content

    @author Name Surname <name.surname@email.org>
    @date 00-00-0000

    A more detailed description of the file (if necessary)
 */

#include <ExampleClass.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
ExampleClass::ExampleClass() :
    M_variableOne (),
    M_variableTwo ()
{

}

ExampleClass::ExampleClass( first_Type&  variableOne,
                            second_Type& variableTwo ) :
    M_variableOne ( variableOne ),
    M_variableTwo ( variableTwo )
{

}

ExampleClass::ExampleClass( const ExampleClass& example ) :
    M_variableOne ( example.M_variableOne ),
    M_variableTwo ( example.M_variableTwo )
{

}

ExampleClass::~ExampleClass()
{

}

// ===================================================
// Operators
// ===================================================
ExampleClass&
ExampleClass::operator=( const ExampleClass& example )
{
    if ( this != &example )
    {
    	M_variableOne = example.M_variableOne;
    	M_variableTwo = example.M_variableTwo;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
ExampleClass::methodOne( first_Type&  inputVariableOne,
                         second_Type& inputVariableTwo )
{
	//Do something
}

void
ExampleClass::methodTwo()
{

}

void
ExampleClass::showMe( std::ostream& output ) const
{
	output << "ExampleClass::showMe()" << std::endl;
	output << "Variable one: " << M_variableOne << std::endl;
	output << "Variable two: " << M_variableTwo << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void
ExampleClass::setVariableOne( const first_Type& variableOne )
{
    M_variableOne = variableOne;
}

// ===================================================
// Get Methods
// ===================================================
const ExampleClass::first_Type&
ExampleClass::variableOne() const
{
    return M_variableOne;
}

// ===================================================
// Private Methods
// ===================================================
void
ExampleClass::privateMethodOne()
{
    //Do something ..
}

} // Namespace LifeV
