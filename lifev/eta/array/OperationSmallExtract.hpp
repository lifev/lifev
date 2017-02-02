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
     @brief This file contains information about the extraction of a row or component within Small* classes.

     @date 08/2012
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef OPERATION_SMALL_EXTRACT_HPP
#define OPERATION_SMALL_EXTRACT_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/MatrixSmall.hpp>

// LifeV namespace.
namespace LifeV
{
//! class OperationSmallExtract  Class containing information about the extraction of a row or component within Small* classes.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class only contains information that can be useful in a templated framework (such as the
  one for assembling the algebraic systems via expressions).

  It cannot be instanciated, neither the generic definition nor the specializations (private constructors
  and destructor only).

  The only information stored in this class is the type of the result, see LifeV::OperationSmallAddition.
*/
template <typename Operand>
class OperationSmallExtract1 // one index specified
{
private:
    //! @name Constructors and destructors
    //@{

    //! No default constructor
    OperationSmallExtract1();

    //! No destructor
    ~OperationSmallExtract1();

    //@}
};

//! \cond


template <UInt Dim1, UInt Dim2>
class OperationSmallExtract1< MatrixSmall<Dim1, Dim2> >
{
public:

    //! @name Public Types
    //@{

    typedef VectorSmall<Dim2> result_Type;

    //@}

private:
    //! @name Constructors and destructors
    //@{

    //! No default constructor
    OperationSmallExtract1();

    //! No destructor
    ~OperationSmallExtract1();

    //@}
};


template <UInt Size>
class OperationSmallExtract1< VectorSmall<Size> >
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
    OperationSmallExtract1();

    //! No destructor
    ~OperationSmallExtract1();

    //@}
};




//! \endcond


template <typename Operand>
class OperationSmallExtract2 // two indexes specified
{
private:
    //! @name Constructors and destructors
    //@{

    //! No default constructor
    OperationSmallExtract2();

    //! No destructor
    ~OperationSmallExtract2();

    //@}
};

//! \cond


template <UInt Dim1, UInt Dim2>
class OperationSmallExtract2< MatrixSmall<Dim1, Dim2> >
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
    OperationSmallExtract2();

    //! No destructor
    ~OperationSmallExtract2();

    //@}
};

//! \endcond

} // namespace LifeV

#endif

