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

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 12 Aug 2011

    A more detailed description of the file (if necessary)
 */

#ifndef TIMEADVANCERYTHMOS_H
#define TIMEADVANCERYTHMOS_H 1

#include <life/lifecore/Life.hpp>

namespace LifeV {

//! TimeAdvanceRythmos - Short description of the class
/*!
    @author Paolo Crosetto
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
class TimeAdvanceRythmos : public TimeAdvance
{
public:

    //! @name Public Types
    //@{

    //@}

    virtual void shiftRight(const feVectorType& solution );
    virtual  void updateRHSFirstDerivative(const Real& timeStep = 1 );
    virtual void updateRHSSecondDerivative(const Real& timeStep = 1 );
    virtual void showMe()  const ;


};

} // Namespace LifeV

#endif /* TIMEADVANCERYTHMOS_H */
