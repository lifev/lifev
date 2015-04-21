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
    @brief File for the definition of the expression used to interpolate FE functions.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_PATCHAREA_HPP
#define EXPRESSION_PATCHAREA_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>

#include <boost/shared_ptr.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionInterpolateValue Class representing an interpolation in an expression.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is meant to be used only in the expression tree, it does not perform
  any computation.

  <b> Template parameters </b>

  <i>MeshType</i>: The type of the mesh

  <i>MapType</i>: The type of the algebraic map (for parallel computations)

  <i>SpaceDim</i>: The ambiant space (for the finite element space)

  <i>FieldDim</i>: The dimension of the field to interpolate (scalar vs vectorial)

  <b> Template requirements </b>

  <i>MeshType</i>: Same as in LifeV::ETFESpace

  <i>MapType</i>: Same as in LifeV::ETFESpace

*/
template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
class ExpressionPatchArea
    : public ExpressionBase<ExpressionPatchArea<MeshType, MapType, SpaceDim, FieldDim>  >
{
public:

    //! @name Public Types
    //@{

    // Base class, used only to make the code cleaner
    typedef ExpressionBase<ExpressionPatchArea<MeshType, MapType, SpaceDim, FieldDim> > base_Type;

    typedef ETFESpace<MeshType, MapType, SpaceDim, FieldDim> ETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace_Type>                ETFESpacePtr_Type;
    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor using the finite element space and the data vector
    ExpressionPatchArea (const ETFESpacePtr_Type feSpace )
        : base_Type(), M_feSpace ( feSpace ) {}

    //! Copy constructor
    ExpressionPatchArea (const ExpressionPatchArea<MeshType, MapType, SpaceDim, FieldDim>& expr)
        : base_Type(), M_feSpace (expr.M_feSpace){}

    //! Destructor
    ~ExpressionPatchArea() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << "patch area" << std::endl;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the finite element space
    const ETFESpacePtr_Type fespace() const
    {
        return M_feSpace;
    }

    // @}

private:

    //! @name Private Methods
    //@{

    //! No default constructor
    ExpressionPatchArea();

    //@}

    // Storage for the finite element space
    ETFESpacePtr_Type M_feSpace;
};

//! Simple function to be used in the construction of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  <b> Template parameters </b>

  <i>MeshType</i>: The type of the mesh stored in the finite element space

  <i>MapType</i>: The type of map used in the finite element space

  <i>SpaceDim</i>: The dimension of the ambient space.

  <i>FieldDim</i>: The dimension of the finite element space (scalar vs vectorial)

  <b> Template requirements </b>

  <i>MeshType</i>: Same as in LifeV::ETFESpace

  <i>MapType</i>: Same as in LifeV::ETFESpace

*/
template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
inline ExpressionPatchArea<MeshType, MapType, SpaceDim, FieldDim>
patchArea (const boost::shared_ptr<ETFESpace<MeshType, MapType, SpaceDim, FieldDim> > fespace)
{
    return ExpressionPatchArea<MeshType, MapType, SpaceDim, FieldDim> ( fespace );
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
