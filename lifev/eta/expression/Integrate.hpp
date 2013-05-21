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
     @brief This file contains the definition of the integrate function.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/util/OpenMPParameters.hpp>

#include <lifev/eta/expression/RequestLoopElement.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/eta/expression/IntegrateMatrixElement.hpp>
#include <lifev/eta/expression/IntegrateVectorElement.hpp>
#include <lifev/eta/expression/IntegrateValueElement.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

/*!
  \namespace ExpressionAssembly

  Namespace for the assembly via expressions

 */
namespace ExpressionAssembly
{

//! Integrate function for matricial expressions
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is an helper function to instantiate the class
  for performing an integration, here to assemble a matrix
  with a loop on the elements.
 */
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression,
            const UInt offsetUp = 0,
            const UInt offsetLeft = 0);
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression,
            const UInt offsetUp,
            const UInt offsetLeft)
{
    return IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>
           (request.mesh(), quadrature, testSpace, solutionSpace, expression, offsetUp, offsetLeft);
}

//! Integrate function for matricial expressions (multi-threaded path)
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is an helper function to instantiate the class
  for performing an integration, here to assemble a matrix
  with a loop on the elements.

  This is an overload of the integrate function for matrices, which
  uses multiple threads to do assembly
 */
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression,
            const OpenMPParameters& ompParams,
            const UInt offsetUp = 0,
            const UInt offsetLeft = 0);
template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression,
            const OpenMPParameters& ompParamsm
            const UInt offsetUp,
            const UInt offsetLeft)
{
    return IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType>
           (request.mesh(), quadrature, testSpace, solutionSpace, expression,
            ompParams, offsetUp, offsetLeft);
}

//! Integrate function for vectorial expressions
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is an helper function to instantiate the class
  for performing an integration, here to assemble a vector
  with a loop on the elements.
 */
template < typename MeshType, typename TestSpaceType, typename ExpressionType>
IntegrateVectorElement<MeshType, TestSpaceType, ExpressionType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const ExpressionType& expression,
            const UInt offset = 0);
template < typename MeshType, typename TestSpaceType, typename ExpressionType>
IntegrateVectorElement<MeshType, TestSpaceType, ExpressionType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<TestSpaceType>& testSpace,
            const ExpressionType& expression,
            const UInt offset)
{
    return IntegrateVectorElement<MeshType, TestSpaceType, ExpressionType>
           (request.mesh(), quadrature, testSpace, expression, offset);
}

//! Integrate function for benchmark expressions
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is an helper function to instantiate the class
  for performing an integration, here to assemble a benchmark
  with a loop on the elements.
 */
template < typename MeshType, typename ExpressionType>
IntegrateValueElement<MeshType, ExpressionType>
integrate ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const ExpressionType& expression)
{
    return IntegrateValueElement<MeshType, ExpressionType>
           (request.mesh(), quadrature, expression);
}


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
