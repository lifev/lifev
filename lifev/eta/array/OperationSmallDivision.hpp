//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains information about the division operation between *Small classes.

     @date 07/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef OPERATION_SMALL_DIVISION_HPP
#define OPERATION_SMALL_DIVISION_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorSmall.hpp>

namespace LifeV
{
//! class OperationSmallDivision  Class containing information about the division operation between *Small classes.

/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    This class only contains information that can be usefull in a templated framework (such as the
    one for assembling the algebraic systems via expressions).

    It cannot be instanciated, neither the generic definition nor the specializations (private constructors
    and destructor only).

    The only information stored in this class is the type of the result, see LifeV::OperationSmallAddition.
*/

template <typename LeftOperand, typename RightOperand>
class OperationSmallDivision
{
private:
    //! @name Constructors and destructors
    //@{

    //! No default constructor
    OperationSmallDivision();

    //! No destructor
    ~OperationSmallDivision();

    //@}
};

//! \cond

template <>
class OperationSmallDivision<Real, Real>
{
public:

    //! @name Public Types
    //@{

    typedef Real result_Type;

    //@}

private:
    //! @name Constructors and destructors
    //@{

    //! No default constructor
    OperationSmallDivision();

    //! No destructor
    ~OperationSmallDivision();

    //@}
};


template <UInt Size>
class OperationSmallDivision< VectorSmall<Size> , Real >
{
public:

    //! @name Public Types
    //@{

    typedef VectorSmall<Size> result_Type;

    //@}

private:
    //! @name Constructors and destructors
    //@{

    //! No default constructor
    OperationSmallDivision();

    //! No destructor
    ~OperationSmallDivision();

    //@}
};

//! \endcond

} // namespace LifeV


#endif

