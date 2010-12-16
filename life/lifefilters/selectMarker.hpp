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
    @brief It contains the standard selector for internal entities

    @author 
    @contributor Nur Aiman Fadel <nur.fadel@mail.polimi.it>
    @maintainer Nur Aiman Fadel <nur.fadel@mail.polimi.it>

    @date

    A more detailed description of the file (if necessary)
 */

#ifndef _SELECTMARKER_HH_
#define _SELECTMARKER_HH_ 1

#include <life/lifemesh/regionMesh3D.hpp>

namespace LifeV
{

//! InternalEntitySelector - Functor class that tells whether an entity flag corresponds to an internal face
/*!
    @author 
    @see 

    This class takes in input an EntityFlag and it can be used <br>
    in order to understand if the input is or not an internal face.
 */

class InternalEntitySelector
{
public:
    //! The default watermark used when standard contructor is adopted.
    static const entityFlag_Type defMarkFlag;

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    InternalEntitySelector();

    //! Short description of the constructor
    /*!
        It is a constructor which requires the EntityFlag
        @param w, the costant EntityFlag which is required in order 
        to create an InternalEntitySelector object.
     */
    InternalEntitySelector(const entityFlag_Type & w);
    //@}


    //! @name Operators
    //@{

    //! The round brackets operator
    /*!
        Operator returning true if the flag corresponds to internal entity.
        If the EntityFlag is greater that the watermark, the associated geometry entity is internal.
        @param test, it is the reference to geometric entity.
        @return true, if the flag corresponds to an internal entity.
     */
    bool operator()(entityFlag_Type const & test) const;
    //@}

private:

    //! The current watermark
    entityFlag_Type M_watermarkFlag;

    //! Unsed operator equivalence
    InternalEntitySelector& operator=( const InternalEntitySelector& example );

}; // Class InternalEntitySelector

} // Namespace LifeV

#endif /* SELECTMARKER_H */
