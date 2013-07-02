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

#ifndef EVALUATE_HPP
#define EVALUATE_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/RequestLoopElement.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/eta/fem/QRAdapterBase.hpp>
#include <lifev/eta/fem/QRAdapterNeverAdapt.hpp>

#include <lifev/eta/expression/EvaluateNodalExpressionVectorElement.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

/*!
  \namespace ExpressionAssembly

  Namespace for the assembly via expressions

 */
namespace ExpressionAssembly
{


//! Integrate function for vectorial expressions
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is an helper function to instantiate the class
  for performing an integration, here to assemble a vector
  with a loop on the elements.
 */
template < typename MeshType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
EvaluateNodalExpressionVectorElement<MeshType, SolutionSpaceType, ExpressionType, QRAdapterType>
evaluateNode ( const RequestLoopElement<MeshType>& request,
            const QRAdapterBase<QRAdapterType>& qrAdapterBase,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression)
{
    return IntegrateVectorElement<MeshType, SolutionSpaceType, ExpressionType, QRAdapterType>
           (request.mesh(), qrAdapterBase.implementation(), solutionSpace, expression);
}

template < typename MeshType, typename SolutionSpaceType, typename ExpressionType>
EvaluateNodalExpressionVectorElement<MeshType, SolutionSpaceType, ExpressionType, QRAdapterNeverAdapt>
evaluateNode ( const RequestLoopElement<MeshType>& request,
            const QuadratureRule& quadrature,
            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
            const ExpressionType& expression)
{
    return IntegrateVectorElement<MeshType, SolutionSpaceType, ExpressionType, QRAdapterNeverAdapt>
           (request.mesh(), QRAdapterNeverAdapt (quadrature), solutionSpace, expression);
}


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
