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
    @brief File where the structures for the element-wise multiplication between expressions are defined.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 08-2012
 */

#ifndef EXPRESSION_DEFORMATIONGRADIENTE_HPP
#define EXPRESSION_DEFORMATIONGRADIENTE_HPP

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/ExpressionBase.hpp>

#include <boost/shared_ptr.hpp>
#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionEmult  Class for representing the transpose operation of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class represents the computation of the deformation gradient matrix operation in the expression tree.

 <b> Template parameters </b>

  <i>MeshType</i>: The type of the mesh stored in the finite element space

  <i>MapType</i>: The type of map used in the finite element space

  <i>SpaceDim</i>: The dimension of the ambient space.

  <i>FieldDim</i>: The dimension of the finite element space (scalar vs vectorial)

  <b> Template requirements </b>

  <i>MeshType</i>: Same as in LifeV::ETFESpace

*/
template<typename MeshType, typename MapType, UInt SpaceDim>
class ExpressionDeformationGradient : public ExpressionBase<ExpressionDeformationGradient<MeshType, MapType, SpaceDim > >
{
public:

    //! @name Public Types
    //@{

    // No real need, just for ease of coding
	typedef ExpressionBase<ExpressionDeformationGradient<MeshType, MapType, SpaceDim> >      base_Type;
    typedef VectorEpetra                                                                     vector_Type;
    typedef boost::shared_ptr<vector_Type>                                                   vectorPtr_Type;
    typedef boost::shared_ptr< ETFESpace<MeshType,MapType,SpaceDim,SpaceDim> >               ETFespacePtr_Type;
    //@}


    //! @name Constructors & Destructor
    //@{

    //! Full constructor
	ExpressionDeformationGradient(const ETFespacePtr_Type fespace, const vectorPtr_Type& vector, const UInt offset) :
        base_Type(),
        M_fespace(fespace),
        M_vector(vector),
        M_offset(offset)
    {}

    //! Copy constructor
	ExpressionDeformationGradient(const ExpressionDeformationGradient<MeshType, MapType, SpaceDim>& expression):
        base_Type(),
        M_fespace(expression.M_fespace),
        M_vector(expression.M_vector),
        M_offset(expression.M_offset)
    {}

    //! Destructor
    ~ExpressionDeformationGradient(){}

    //@}


    //! @name Methods
    //@{

    //! Display method
	static void display(std::ostream& out= std::cout)
    { out << " Deformation Gradient "; }

    //@}


    //! @name Get Methods
    //@{
	const vector_Type& vector() const {return *M_vector;}

	const UInt offset() const {return M_offset;}

    ETFespacePtr_Type fespace() const {return M_fespace;}
    //@}

private:

    //! @name Private Methods
    //@{

    //! No default constructor
	ExpressionDeformationGradient();

    //@}

    // Storage for the finite element space
	ETFespacePtr_Type M_fespace;

    // The vector to compute the deformation gradient
	vectorPtr_Type    M_vector;
    UInt              M_offset;
};


// transpose  The generic function for the transpose of an expression.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  Function used in construction of the expression tree. To avoid shadowing
  other functions, it uses the ExpressionBase type to distinguish expressions
  from other types.

  <b> Template parameters </b>

  <i>ExpressionType</i>: The expression to be transposed.

  <b> Template requirements </b>

  <i>ExpressionType</i>: Same as in LifeV::ExpressionDeformationGradient

*/
template<typename MeshType, typename MapType, UInt SpaceDim>
ExpressionDeformationGradient<MeshType, MapType, SpaceDim>
F(const boost::shared_ptr< ETFESpace<MeshType,MapType,SpaceDim,SpaceDim> > fespace,
  const boost::shared_ptr<VectorEpetra>& vectorPtr,
  UInt offset)
{
	return ExpressionDeformationGradient<MeshType, MapType, SpaceDim>(fespace, vectorPtr, offset);
};

} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
