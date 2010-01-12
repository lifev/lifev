//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-(>>>YEAR<<<) EPFL, Politecnico di Milano, INRIA

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
    @file
    @brief A short description of the file content

    @author (>>>USER_NAME<<<) <(>>>AUTHOR<<<)>
    @date (>>>DATE<<<)

    A more detailed description of the file (if necessary)
 */

#include <(>>>FILE_SANS<<<).hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
(>>>FILE_SANS<<<)::(>>>FILE_SANS<<<)() :
    M_variableOne (),
    M_variableTwo ()
{

}

(>>>FILE_SANS<<<)::(>>>FILE_SANS<<<)( first_Type&  variableOne,
                              second_Type& variableTwo ) :
    M_variableOne ( variableOne ),
    M_variableTwo ( variableTwo )
{

}

(>>>FILE_SANS<<<)::(>>>FILE_SANS<<<)( const TemplateClass& T ) :
    M_variableOne ( T.M_variableOne ),
    M_variableTwo ( T.M_variableTwo )
{

}

(>>>FILE_SANS<<<)::~(>>>FILE_SANS<<<)()
{

}

// ===================================================
// Operators
// ===================================================
(>>>FILE_SANS<<<)&
(>>>FILE_SANS<<<)::operator=( const (>>>FILE_SANS<<<)& T )
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
(>>>FILE_SANS<<<)::methodOne( first_Type&  inputVariableOne,
                          second_Type& inputVariableTwo )
{
	//Do something
}

void
(>>>FILE_SANS<<<)::methodTwo()
{

}

void
(>>>FILE_SANS<<<)::showMe( std::ostream& output ) const
{
	output << "TemplateClass::showMe()" << std::endl;
	output << "Variable one: " << M_variableOne << std::endl;
	output << "Variable two: " << M_variableTwo << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void
(>>>FILE_SANS<<<)::setVariableOne( const first_Type& variableOne )
{
    M_variableOne = variableOne;
}

// ===================================================
// Get Methods
// ===================================================
const (>>>FILE_SANS<<<)::first_Type&
(>>>FILE_SANS<<<)::variableOne() const
{
    return M_variableOne;
}

// ===================================================
// Private Methods
// ===================================================
void
(>>>FILE_SANS<<<)::privateMethodOne()
{
    //Do something ..
}

} // Namespace LifeV
