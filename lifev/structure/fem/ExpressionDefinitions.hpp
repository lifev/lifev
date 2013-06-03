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
    @brief File containing the procedures for the local assembly of the differential operators

    @author Gianmarco Mengaldo <gianmarco.mengaldo@gmail.com>
    @author Paolo Tricerri <gianmarco.mengaldo@gmail.com>
    @mantainer Paolo Tricerri <paolo.tricerri@epfl.ch>

    All the methods are described in the report StructuralSolver framework in LifeV: Description and Usage
 */


#ifndef _EXPRESSIONDEFINITIONS_H_
#define _EXPRESSIONDEFINITIONS_H_

#include <string>
#include <iostream>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV
{

typedef RegionMesh<LinearTetra>                       MeshType;
typedef ETFESpace<MeshType, MapEpetra, 3, 3 >         ETFESpace_Type;
typedef VectorEpetra                                  vector_Type;
typedef MatrixSmall<3,3>                              matrixSmall_Type;

/*! /namespace ExpressionDefinitions

  This namespace is specially designed to contain the elementary
  operations (corresponding to differential operators) that build
  the local contributions to be used in the assembly procedures.

*/
namespace ExpressionDefinitions
{

using namespace ExpressionAssembly;

//! @name Public typedefs
//@{

typedef ExpressionAddition<
    ExpressionInterpolateGradient<MeshType, MapEpetra, 3, 3>, ExpressionMatrix<3,3> >  deformationGradient_Type;

typedef ExpressionDeterminant<
    ExpressionAddition< ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > determinantTensorF_Type;

typedef ExpressionProduct<
    ExpressionTranspose<
        ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
    ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> >
    > rightCauchyGreenTensor_Type;

typedef ExpressionMinusTransposed<
    ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > minusTransposedTensor_Type;

typedef ExpressionTrace<
    ExpressionProduct<ExpressionTranspose<ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
                          ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >
    > traceTensor_Type;

typedef ExpressionDot<
  ExpressionProduct<
    ExpressionTranspose<
      ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
    ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
  ExpressionProduct<
    ExpressionTranspose<
      ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
    ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > > traceSquaredTensor_Type;

typedef ExpressionPower<
  ExpressionDeterminant<
    ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > > powerExpression_Type;

typedef ExpressionProduct<
  ExpressionPower<
    ExpressionDeterminant< ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > >,
  ExpressionTrace<
    ExpressionProduct<
      ExpressionTranspose<ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
      ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > > > isochoricTrace_Type;


  // Typedefs for anisotropic laws
#ifdef ENABLE_ANISOTROPIC_LAW
  typedef  ExpressionInterpolateValue<MeshType, MapEpetra, 3, 3>   interpolatedValue_Type;

  typedef  ExpressionOuterProduct<
    ExpressionInterpolateValue<MeshType, MapEpetra, 3, 3>,
    ExpressionInterpolateValue<MeshType, MapEpetra, 3, 3> >        outerProduct_Type;

  typedef  ExpressionDot<
    ExpressionProduct<
      ExpressionTranspose<
	ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
      ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
    ExpressionOuterProduct<
      ExpressionInterpolateValue<MeshType, MapEpetra, 3, 3 >, ExpressionInterpolateValue<MeshType, MapEpetra, 3, 3 > > > stretch_Type;

  typedef  ExpressionProduct<
    ExpressionPower<
      ExpressionDeterminant<
	ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > >,

    ExpressionDot<
      ExpressionProduct<
	ExpressionTranspose<
	  ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
	ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
      ExpressionOuterProduct<
	ExpressionInterpolateValue<MeshType, MapEpetra, 3, 3 >, ExpressionInterpolateValue<MeshType, MapEpetra, 3, 3 > > > > isochoricStretch_Type;

  // This typedef describes \kappa * trCBar
  typedef ExpressionProduct< ExpressionScalar, isochoricTrace_Type > distributedIsochoricTrace_Type;

  // This typedef describes ( 1 - 3 \kappa) * IVithBar 
  typedef ExpressionProduct< ExpressionScalar, isochoricStretch_Type > distributedIsochoricStretch_Type;

  // This typedef describes \kappa * trCBar + ( 1 - 3 \kappa) * IVithBar 
  typedef ExpressionAddition< 
    distributedIsochoricTrace_Type, 
    distributedIsochoricStretch_Type > distributedInvariants_Type;

  // This typedef describes \kappa * trCBar + ( 1 - 3 \kappa) * IVithBar - 1.0
  typedef ExpressionAddition<
    distributedInvariants_Type,
    ExpressionScalar>                                              distributedStretch_Type;

#endif

//@}

deformationGradient_Type deformationGradient( const boost::shared_ptr< ETFESpace_Type > dispETFESpace,
                                              const vector_Type& disp, UInt offset, const matrixSmall_Type identity)
{
    return deformationGradient_Type( grad( dispETFESpace,  disp, offset), value(identity) );
}

determinantTensorF_Type determinantF( const deformationGradient_Type F )
{
    return determinantTensorF_Type( F );
}

rightCauchyGreenTensor_Type tensorC( const ExpressionTranspose<deformationGradient_Type> tF, const deformationGradient_Type F )
{
    return rightCauchyGreenTensor_Type( tF, F );
}

minusTransposedTensor_Type minusT( const deformationGradient_Type F )
{
    return minusTransposedTensor_Type( F );
}

traceTensor_Type traceTensor( const rightCauchyGreenTensor_Type C )
{
    return traceTensor_Type( C );
}

traceSquaredTensor_Type traceSquared( const rightCauchyGreenTensor_Type C )
{
  return traceSquaredTensor_Type( C , C );
}

powerExpression_Type powerExpression( const determinantTensorF_Type J, const Real exponent )
{
  return powerExpression_Type( J , exponent );
}

isochoricTrace_Type isochoricTrace( const powerExpression_Type Jel, const traceTensor_Type I )
{
  return isochoricTrace_Type( Jel , I );
}



// Constructors for anisotropic laws
#ifdef ENABLE_ANISOTROPIC_LAW
interpolatedValue_Type interpolateFiber( const boost::shared_ptr< ETFESpace_Type > dispETFESpace,
					 const vector_Type& fiberVector)
{
  return interpolatedValue_Type( dispETFESpace, fiberVector ) ;
}

outerProduct_Type fiberTensor( const interpolatedValue_Type ithFiber )
{
  return outerProduct_Type( ithFiber, ithFiber );
}

stretch_Type fiberStretch( const rightCauchyGreenTensor_Type C, const outerProduct_Type M)
{
  return stretch_Type( C, M );
}

isochoricStretch_Type isochoricFourthInvariant( const powerExpression_Type Jel, const stretch_Type I_4ith)
{
  return isochoricStretch_Type( Jel, I_4ith );
}

distributedStretch_Type distributedStretch( const isochoricTrace_Type trCBar, const isochoricStretch_Type I_4ith, const Real kappa )
{
  distributedIsochoricTrace_Type dIC_bar( value(kappa), trCBar ) ;
  distributedIsochoricStretch_Type dI4_bar( value( 1.0 - 3.0 * kappa), I_4ith ) ;

  distributedInvariants_Type dInvariants( dIC_bar, dI4_bar );

  return distributedStretch_Type( dInvariants, value( -1.0 ) );
}

#endif


} //! End namespace ExpressionDefinitions

} //! End namespace LifeV
#endif
