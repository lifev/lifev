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
    @brief File containing the definition of the scalar expression.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_SCALAR_HPP
#define EXPRESSION_SCALAR_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionScalar  Class representing a scalar constant in an expression.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
*/
class ExpressionScalar : public ExpressionBase<ExpressionScalar>
{
public:

    //! @name Public Types
    //@{

    typedef ExpressionBase<ExpressionScalar> base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor using the value of the scalar
    ExpressionScalar (const Real& myValue);

    //! Copy constructor
    ExpressionScalar (const ExpressionScalar& expr);

    //! Destructor
    ~ExpressionScalar();

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout);

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the value of the scalar
    const Real& value() const;

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    ExpressionScalar();

    //@}

    Real M_value;
};

//! Simple function to be used in the construction of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
*/
inline ExpressionScalar
value (const Real& myValue)
{
    return ExpressionScalar (myValue);
}


//! class ExpressionScalar  Class representing a scalar constant in an expression.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
*/

template<typename VectorType>
class ExpressionExtractScalar : public ExpressionBase<ExpressionExtractScalar<VectorType> >
{
public:

    //! @name Public Types
    //@{

	typedef ExpressionBase<ExpressionExtractScalar<VectorType> > base_Type;
    typedef boost::shared_ptr<VectorType>                        containerPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor using the value of the scalar
	ExpressionExtractScalar(const VectorType& myVector)
    : base_Type(), M_vector(new VectorType(myVector) ) {}

    //! Copy constructor
	ExpressionExtractScalar(const ExpressionExtractScalar<VectorType>& expr)
    : base_Type(), M_vector( expr.vector() ){}

    //! Destructor
    ~ExpressionExtractScalar(){}

    //@}


    //! @name Methods
    //@{

    //! Display method
	static void display(std::ostream& out= std::cout);

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the value of the scalar
	const containerPtr_Type vector() const { return M_vector; } ;

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    ExpressionExtractScalar();

    //@}

	containerPtr_Type M_vector;
};

//! Simple function to be used in the construction of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
*/
template<typename VectorType>
inline ExpressionExtractScalar<VectorType>
parameter(const VectorType& myVector)
{
	return ExpressionExtractScalar<VectorType>( myVector );
}

} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
