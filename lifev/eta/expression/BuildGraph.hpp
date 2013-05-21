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
     @brief This file contains the definition of the buildGraph function.

     This function is used to precompute the graph of a finite element
     matrix, allowing the matrix to be build in closed, optimized form,
     which makes the assembly procedure more efficient

     @date 03/2013
     @author Radu Popescu <radu.popescu@epfl.ch>
 */

#ifndef BUILD_GRAPH_HPP
#define BUILD_GRAPH_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/RequestLoopElement.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/eta/expression/GraphElement.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

/*!
  \namespace ExpressionAssembly

  Namespace for the assembly via expressions

 */
namespace ExpressionAssembly
{

//! Function to precompute the matrix graph
/*!
  @author Radu Popescu <radu.popescu@epfl.ch>

  This is a helper function to precompute the Crs graph used to build
  a FECrsMatrix in closed optimized form
 */
template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType >
GraphElement < MeshType,
             TestSpaceType,
             SolutionSpaceType,
             ExpressionType >
             buildGraph ( const RequestLoopElement<MeshType>& request,
                          const QuadratureRule& quadrature,
                          const boost::shared_ptr<TestSpaceType>& testSpace,
                          const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
                          const ExpressionType& expression,
                          const UInt offsetUp = 0,
                          const UInt offsetLeft = 0);
template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType >
GraphElement < MeshType,
             TestSpaceType,
             SolutionSpaceType,
             ExpressionType >
             buildGraph ( const RequestLoopElement<MeshType>& request,
                          const QuadratureRule& quadrature,
                          const boost::shared_ptr<TestSpaceType>& testSpace,
                          const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
                          const ExpressionType& expression,
                          const UInt offsetUp,
                          const UInt offsetLeft)
{
    return GraphElement < MeshType,
           TestSpaceType,
           SolutionSpaceType,
           ExpressionType >
           (request.mesh(), quadrature, testSpace, solutionSpace, expression,
            offsetUp, offsetLeft );
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif // BUILD_GRAPH_HPP
