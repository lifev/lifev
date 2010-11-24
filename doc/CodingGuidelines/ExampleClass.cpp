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

#include <TemplateClass.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
TemplateClass::TemplateClass() :
    M_variableOne (),
    M_variableTwo ()
{

}

TemplateClass::TemplateClass( first_Type&  variableOne,
                              second_Type& variableTwo ) :
    M_variableOne ( variableOne ),
    M_variableTwo ( variableTwo )
{

}

TemplateClass::TemplateClass( const TemplateClass& T ) :
    M_variableOne ( T.M_variableOne ),
    M_variableTwo ( T.M_variableTwo )
{

}

TemplateClass::~TemplateClass()
{

}

// ===================================================
// Operators
// ===================================================
TemplateClass&
TemplateClass::operator=( const TemplateClass& T )
{
    if ( this != &T )
    {
    	M_variableOne = T.M_variableOne;
    	M_variableTwo = T.M_variableTwo;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
TemplateClass::methodOne( first_Type&  inputVariableOne,
                          second_Type& inputVariableTwo )
{
	//Do something
}

void
TemplateClass::methodTwo()
{

}

void
TemplateClass::showMe( std::ostream& output ) const
{
	output << "TemplateClass::showMe()" << std::endl;
	output << "Variable one: " << M_variableOne << std::endl;
	output << "Variable two: " << M_variableTwo << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void
TemplateClass::setVariableOne( const first_Type& variableOne )
{
    M_variableOne = variableOne;
}

// ===================================================
// Get Methods
// ===================================================
const TemplateClass::first_Type&
TemplateClass::variableOne() const
{
    return M_variableOne;
}

// ===================================================
// Private Methods
// ===================================================
void
TemplateClass::privateMethodOne()
{
    //Do something ..
}

} // Namespace LifeV
