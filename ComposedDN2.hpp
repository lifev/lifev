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

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 08 Jun 2010

    A more detailed description of the file (if necessary)
 */

#ifndef COMPOSEDDN2_H
#define COMPOSEDDN2_H 1

#include <lifemc/lifesolver/ComposedDN.hpp>

namespace LifeV {

//! ComposedDN - Short description of the class
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
class ComposedDN2 : public ComposedDN
{
public:
    typedef ComposedDN super;

    ComposedDN2():
        super(11)
    {
    }

    void coupler(map_shared_ptrtype map,
                 const std::map<ID, ID>& locDofMap,
                 const vector_ptrtype numerationInterface,
                 const Real& timeStep);

    virtual void replace_matrix( const matrix_ptrtype& oper, UInt position);
    virtual void replace_bch(bchandler_ptrtype& oper, UInt position){M_bch[1-position]=oper;}
    void blockAssembling();

private:

    void replace_precs( matrix_ptrtype& oper, UInt position);
    //    static bool                                      reg;

};

} // Namespace LifeV

#endif /* COMPOSEDDN2_H */
