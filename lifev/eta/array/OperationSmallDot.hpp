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
     @brief This file contains information about the dot product operation between Small* classes.

     @date 07/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef OPERATION_SMALL_DOT_HPP
#define OPERATION_SMALL_DOT_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorSmall.hpp>

// LifeV namespace.
namespace LifeV
{
//! class OperationSmallDot  Class containing information about the dot product operation between Small* classes.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class only contains information that can be usefull in a templated framework (such as the
  one for assembling the algebraic systems via expressions).

  It cannot be instanciated, neither the generic definition nor the specializations (private constructors
  and destructor only).

  The only information stored in this class is the type of the result, see LifeV::OperationSmallAddition.
*/

template <typename LeftOperand, typename RightOperand>
class OperationSmallDot
{
private:
    //! @name Constructors and destructors
    //@{

    //! No default constructor
    OperationSmallDot();

    //! No destructor
    ~OperationSmallDot();

    //@}
};

//! \cond

template <UInt Size>
class OperationSmallDot< VectorSmall<Size> , VectorSmall<Size> >
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
    OperationSmallDot();

    //! No destructor
    ~OperationSmallDot();

    //@}
};

template <UInt Size1, UInt Size2>
class OperationSmallDot< MatrixSmall<Size1, Size2> , MatrixSmall<Size1, Size2> >
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
    OperationSmallDot();

    //! No destructor
    ~OperationSmallDot();

    //@}
};

// specialization for matrix * vector
template <UInt Dim1, UInt Dim2>
class OperationSmallDot< MatrixSmall<Dim1, Dim2> , VectorSmall<Dim2> >
{
public:

    //! @name Public Types
    //@{

    typedef VectorSmall<Dim1> result_Type;

    //@}

private:
    //! @name Constructors and destructors
    //@{

    //! No default constructor
    OperationSmallDot();

    //! No destructor
    ~OperationSmallDot();

    //@}
};

// specialization for matrix * vector
template <UInt Dim1, UInt Dim2>
class OperationSmallDot< VectorSmall<Dim2>, MatrixSmall<Dim1, Dim2> >
{
public:

    //! @name Public Types
    //@{

    typedef VectorSmall<Dim1> result_Type;

    //@}

private:
    //! @name Constructors and destructors
    //@{

    //! No default constructor
    OperationSmallDot();

    //! No destructor
    ~OperationSmallDot();

    //@}
};

//! \endcond

} // namespace LifeV

#endif

