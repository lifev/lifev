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
    @contributor Name Surname <name.surname@email.org>
    @maintainer Name Surname <name.surname@email.org>

    @date 00-00-0000

    A more detailed description of the file (if necessary)
 */

#ifndef EXAMPLECLASS_H
#define EXAMPLECLASS_H 1

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

//! ExampleClass - Short description of the class
/*!
    @author Name Surname
    @see Reference to papers (if available)

    Here write a long and detailed description of the class.

    For this purpose you can use a lot of standard HTML code.
    Here there is a list with some useful examples:

    For bold text use: <b>BOLD TEXT</b>

    For empatyze a word type @e followed by the word

    For verbatim a word type @c followed by the word

    For vertical space (empty lines) use: <br>

    For creating list type:
    <ol>
        <li> First element of the enumerated list
        <ul>
             <li> First element of the dotted sublist.
             <li> Second element of the dotted sublist
        </ul>
        <li> Second element of the enumerated list
        <li> Third element of the enumerated list
    </ol>

    For writing a warning type: @warning followed by the description
    of the warning

    It is possible to use a lot of other standard HTML commands.
    Visit http://www.stack.nl/~dimitri/doxygen/htmlcmds.html for
    a detailed list.

    For any other kind of information visit www.doxygen.org.
 */
class ExampleClass
{
public:

    //! @name Public Types
    //@{

    /*! @enum listOfTemplatesOptions
        Description of the purpose of the enumerator list.
     */
    enum listOfTemplatesOptions
    {
        options1, /*!< This options means ... */
        options2, /*!< This options means ... */
        options3  /*!< This options means ... */
    };

    typedef int                                first_Type;
    typedef double                             second_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    ExampleClass();

    //! Short description of the constructor
    /*!
        Add more details about the constructor.
        NOTE: short description is automatically added before this part.
        @param VariableOne Description of the first variable
        @param VariableTwo Description of the second variable
     */
    ExampleClass ( first_Type& VariableOne, second_Type& VariableTwo );

    //! Copy constructor
    /*!
        Add more details about the copy constructor.
        NOTE: short description is automatically added before this part.
        @param example ExampleClass
     */
    ExampleClass ( const ExampleClass& example );

    //! Destructor
    ~ExampleClass();

    //@}


    //! @name Operators
    //@{

    //! The equivalence operator
    /*!
        Add more details about the method.
        NOTE: short description is automatically added before this part.
        @param example ExampleClass
        @return Reference to a new ExampleClass with the same
                content of ExampleClass example
     */
    ExampleClass& operator= ( const ExampleClass& example );

    //@}


    //! @name Methods
    //@{

    //! Short description of this method
    /*!
        Add more details about the method.
        NOTE: short description is automatically added before this part.
        @param inputVariableOne Description of the first input variable
        @param inputVariableTwo Description of the second input variable
     */
    void methodOne ( first_Type& inputVariableOne, second_Type& inputVariableTwo );

    //! Short description of this method
    /*!
        Add more details about the method.
        NOTE: short description is automatically added before this part.
     */
    void methodTwo();

    //! Display general information about the content of the class
    /*!
        List of things displayed in the class
        @param output specify the output format (std::cout by default)
     */
    void showMe ( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Short description of this set method
    /*!
        Add more details about the method.
        NOTE: short description is automatically added before this part.
        @param variableOne Description of the input variable
     */
    void setVariableOne ( const first_Type& variableOne );

    //@}


    //! @name Get Methods
    //@{

    //! Short description of this get method
    /*!
        Add more details about the method.
        NOTE: short description is automatically added before this part.
        @return Description of the output variable
     */
    const first_Type& variableOne() const;

    //@}

private:

    //! @name Private Methods
    //@{

    //! Short description of this method
    /*!
        Add more details about the method.
        NOTE: short description is automatically added before this part.
     */
    void privateMethodOne();

    //@}

    first_Type                       M_variableOne;
    second_Type                      M_variableTwo;
};

} // Namespace LifeV

#endif /* EXAMPLECLASS_H */
