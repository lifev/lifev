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
    @file
    @brief A short description of the file content

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 25 Nov 2010

    A more detailed description of the file (if necessary)
 */

#ifndef QUADRATURERULEPROVIDER_H
#define QUADRATURERULEPROVIDER_H 1

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>

namespace LifeV
{

//! QuadratureRuleProvider - Short description of the class
/*!
    @author Samuel Quinodoz
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
class QuadratureRuleProvider
{
public:

    //! @name Public Types
    //@{

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Destructor
    virtual ~QuadratureRuleProvider();

    //@}



    //! @name Methods
    //@{

    //! Provide a quadrature rule with the given exactness (abort if not available)
    /*!
      Given a shape, this method will try to return a quadrature rule that has the
      given exactness. If such a quadrature rule is not defined, the program will
      abort.
     */
    static const QuadratureRule& provideExactness (const ReferenceShapes& shape, const UInt& exactness);

    //! Provide the quadrature rule with the highest exactness available.
    static const QuadratureRule& provideMaximal (const ReferenceShapes& shape);

    //! Provide a quadrature with the given exactness (or the maximal one if not available)
    /*
      Given a shape, this method tries to return a quadrature with the given exactness.
      If such a quadrature is not defined, the quadrature with the highest exactness for
      the shape is returned.
     */
    static const QuadratureRule& provideExactnessMax (const ReferenceShapes& shape, const UInt& exactness);

    //@}



    //! @name Operators
    //@{

    //@}


    //! @name Set Methods
    //@{

    //@}


    //! @name Get Methods
    //@{

    //@}

private:

    //! @name Private Methods
    //@{

    //! Empty Constructor
    QuadratureRuleProvider();

    //! Copy Constructor
    QuadratureRuleProvider ( const QuadratureRuleProvider& T );

    //! Method for the differentShapes

    static const QuadratureRule& provideExactnessTetra (const UInt& exactness);
    static const QuadratureRule& provideExactnessPrism (const UInt& exactness);
    static const QuadratureRule& provideExactnessHexa (const UInt& exactness);
    static const QuadratureRule& provideExactnessQuad (const UInt& exactness);
    static const QuadratureRule& provideExactnessTriangle (const UInt& exactness);
    static const QuadratureRule& provideExactnessLine (const UInt& exactness);
    static const QuadratureRule& provideExactnessPoint (const UInt& exactness);

    //@}

    static const UInt S_maxExactnessTetra;
    static const UInt S_maxExactnessPrism;
    static const UInt S_maxExactnessHexa;
    static const UInt S_maxExactnessQuad;
    static const UInt S_maxExactnessTriangle;
    static const UInt S_maxExactnessLine;
    // No static for point, it is always exact!

};

} // Namespace LifeV

#endif /* QUADRATURERULEPROVIDER_H */
