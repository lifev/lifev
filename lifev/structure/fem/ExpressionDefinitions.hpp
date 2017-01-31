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

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV
{

typedef RegionMesh<LinearTetra>                       MeshType;
typedef ETFESpace<MeshType, MapEpetra, 3, 3 >         ETFESpace_Type;
typedef ETFESpace<MeshType, MapEpetra, 3, 1 >         scalarETFESpace_Type;
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

// Definition of F = \grad(displacement) + I
typedef ExpressionAddition<
    ExpressionInterpolateGradient<MeshType, MapEpetra, 3, 3>, ExpressionMatrix<3,3> >  deformationGradient_Type;

// Definition of J = det(F)
typedef ExpressionDeterminant<
    ExpressionAddition< ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > determinantTensorF_Type;

// Definition of the right Cauchy-Green tensor C = F^{T} * F
typedef ExpressionProduct<
    ExpressionTranspose<
        ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
    ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> >
    > rightCauchyGreenTensor_Type;

// Definition of the tensor F^{-T}
typedef ExpressionMinusTransposed<
    ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > minusTransposedTensor_Type;

// Used later specially in the multi-mechanism
// Definition of the tensor F^{-1}
typedef ExpressionInverse<
    ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > inverseTensor_Type;


// Definition of the trace of the tensor C
typedef ExpressionTrace<
    ExpressionProduct<ExpressionTranspose<ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
                          ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >
    > traceTensor_Type;

// Definition of the trace of the tensor C^2 = C:C = tr(C^2)
typedef ExpressionDot<
  ExpressionProduct<
    ExpressionTranspose<
      ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
    ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
  ExpressionProduct<
    ExpressionTranspose<
      ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
    ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > > traceSquaredTensor_Type;

// Definition of the power of J ( specifically J^(-2.0/3.0) )
typedef ExpressionPower<
  ExpressionDeterminant<
    ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > > powerExpression_Type;

// Definition of the power of J ( specifically J^(-2.0/3.0) )
typedef ExpressionIsochoricChangeOfVariable< determinantTensorF_Type> isochoricChangeOfVariable_Type;

// Definition of the isochoric trace \bar{I_C} = J^( -2.0/3.0)*tr(C))
typedef ExpressionProduct<
  ExpressionPower<
    ExpressionDeterminant< ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > >,
  ExpressionTrace<
    ExpressionProduct<
      ExpressionTranspose<ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > >,
      ExpressionAddition<ExpressionInterpolateGradient<MeshType, MapEpetra,3,3>, ExpressionMatrix<3,3> > > > > isochoricTrace_Type;


  // Typedefs for anisotropic laws
  typedef  ExpressionInterpolateValue<MeshType, MapEpetra, 3, 3>   interpolatedValue_Type;

  typedef  ExpressionInterpolateValue<MeshType, MapEpetra, 3, 1>   interpolatedScalarValue_Type;

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

//@}

 deformationGradient_Type deformationGradient( const std::shared_ptr< ETFESpace_Type > dispETFESpace,
					       const vector_Type& disp, UInt offset, const matrixSmall_Type identity);

  determinantTensorF_Type determinantF( const deformationGradient_Type F );

  rightCauchyGreenTensor_Type tensorC( const ExpressionTranspose<deformationGradient_Type> tF, const deformationGradient_Type F );

  minusTransposedTensor_Type minusT( const deformationGradient_Type F );

  inverseTensor_Type inv( const deformationGradient_Type F );

  traceTensor_Type traceTensor( const rightCauchyGreenTensor_Type C );

  traceSquaredTensor_Type traceSquared( const rightCauchyGreenTensor_Type C );

  powerExpression_Type powerExpression( const determinantTensorF_Type J, const Real exponent );

  isochoricChangeOfVariable_Type isochoricDeterminant( const determinantTensorF_Type J );

  isochoricTrace_Type isochoricTrace( const powerExpression_Type Jel, const traceTensor_Type I );

// Constructors for anisotropic laws
  interpolatedValue_Type interpolateFiber( const std::shared_ptr< ETFESpace_Type > dispETFESpace,
					  const vector_Type& fiberVector);

  interpolatedValue_Type interpolateValue( const std::shared_ptr< ETFESpace_Type > dispETFESpace,
					   const vector_Type& valueVector);

  interpolatedScalarValue_Type interpolateScalarValue( const std::shared_ptr< scalarETFESpace_Type > dispETFESpace,
                                                       const vector_Type& valueVector);

  outerProduct_Type fiberTensor( const interpolatedValue_Type ithFiber );

  stretch_Type fiberStretch( const rightCauchyGreenTensor_Type C, const outerProduct_Type M);

  isochoricStretch_Type isochoricFourthInvariant( const powerExpression_Type Jel, const stretch_Type I_4ith);
} //! End namespace ExpressionDefinitions


//! The namespace ExpressionDistributedModel is specific for the Dstributed Holzapfel model
//! the definitions have been inserted here in order to avoid huge declarations of expressions
//! in the header file of the model

//! The goal of the current namespace is just to define the expressions (final and intermediate)
//! that are needed for the distributed model making use of the previous namespaces already
//! defined.

namespace ExpressionDistributedModel
{
using namespace ExpressionAssembly;

//! @name Public typedefs
//@{

// Definition of the expression which represents
// the derivative with respect to F of the distributed
// stretch of the fibers.

// This typedef describes \kappa * trCBar
typedef ExpressionProduct< ExpressionScalar, ExpressionDefinitions::isochoricTrace_Type > distributedIsochoricTrace_Type;

// This typedef describes ( 1 - 3 \kappa) * IVithBar
typedef ExpressionProduct< ExpressionScalar, ExpressionDefinitions::isochoricStretch_Type > distributedIsochoricStretch_Type;

// This typedef describes \kappa * trCBar + ( 1 - 3 \kappa) * IVithBar
typedef ExpressionAddition<
    distributedIsochoricTrace_Type,
    distributedIsochoricStretch_Type > distributedInvariants_Type;

// This typedef describes \kappa * trCBar + ( 1 - 3 \kappa) * IVithBar - 1.0
typedef ExpressionAddition<
    distributedInvariants_Type,
    ExpressionScalar>                                              distributedStretch_Type;

// Term that represents F^-T : dF
typedef ExpressionDot<
    ExpressionDefinitions::minusTransposedTensor_Type, ExpressionDphiJ > minusTFscalarDF_distrType;

// Term that represents F : dF
typedef ExpressionDot<
    ExpressionDefinitions::deformationGradient_Type, ExpressionDphiJ > FscalarDF_distrType;

// Term that represents dF^T F : M
typedef ExpressionDot<
    ExpressionProduct< ExpressionTranspose<ExpressionDphiJ>, ExpressionDefinitions::deformationGradient_Type >,
    ExpressionDefinitions::outerProduct_Type  >                        dFTtimesFscalarM_distrType;

// Term that represents F^T dF : M
typedef ExpressionDot<
    ExpressionProduct< ExpressionTranspose<ExpressionDefinitions::deformationGradient_Type>, ExpressionDphiJ >,
    ExpressionDefinitions::outerProduct_Type  >                        FTtimesDFscalarM_distrType;

/*
  Term that represents the first derivative of the isochoric trace with respect to F
  In formula:

  D_F( \bar(I_C) ) : dF = D_F( J^(-2.0/3.0) * I_C ) : dF =
  ( -2.0/3.0 ) * \bar(I_C) F^{-T}:dF + 2 J^(-2.0/3.0) F:dF
*/

typedef ExpressionProduct< ExpressionScalar, ExpressionDefinitions::traceTensor_Type>           scaledTrace_Type;
typedef ExpressionProduct< ExpressionScalar, ExpressionDefinitions::isochoricTrace_Type>        scaledIsochoricTrace_Type;
typedef ExpressionProduct< ExpressionScalar, ExpressionDefinitions::powerExpression_Type>       scaledDeterminant_Type;

typedef ExpressionProduct< ExpressionScalar, ExpressionDefinitions::stretch_Type>               scaledFourthInvariant_Type;
typedef ExpressionProduct< ExpressionScalar, ExpressionDefinitions::isochoricStretch_Type>      scaledIsochoricFourthInvariant_Type;

typedef ExpressionAddition<
    ExpressionProduct< scaledIsochoricTrace_Type, minusTFscalarDF_distrType>,
    ExpressionProduct< scaledDeterminant_Type,    FscalarDF_distrType>
    > linearizationFisochoricTrace_Type;

/*
  Term that represents the first derivative of the isochoric fourth invariant
  In formula:

  D_F( \bar(I4^i) ) : dF = D_F( J^(-2.0/3.0) * I4^i ) : dF =
  ( -2.0/3.0 ) * \bar(I4^i) F^{-T}:dF + J^(-2.0/3.0) ( dF^T F + F^T dF ): M
*/
typedef ExpressionAddition<
    ExpressionProduct< scaledIsochoricFourthInvariant_Type, minusTFscalarDF_distrType>,
    ExpressionProduct< ExpressionDefinitions::powerExpression_Type,
                       ExpressionAddition<dFTtimesFscalarM_distrType, FTtimesDFscalarM_distrType > >
    > linearizationFisochoricFourthInvariant_Type;

/*
  Term that represents the sum of the two linearizations to have the derivative
  of the distributed stretch
*/

typedef ExpressionAddition<
    ExpressionProduct< ExpressionScalar, linearizationFisochoricTrace_Type > ,
    ExpressionProduct< ExpressionScalar, linearizationFisochoricFourthInvariant_Type >
    > linearizationDistributedStretch_Type;

//====================================================

/*
  Definition of the expression that is the tensorial part
  of the first Piola-Kirchhoff tensor for this model

  In the derivation of the Piola-Kirchhoff tensor,
  the tensorial part is of the form

  \kappa F - ( 1/3 )\kappa I_C F^{-T} +
  + ( 1 - 3\kappa ) FM - ( (1-3\kappa) / 3 ) * I_4 * F^{-T}
*/
typedef ExpressionProduct< ExpressionScalar, ExpressionDefinitions::deformationGradient_Type>    scaledTensorF_Type;
typedef ExpressionProduct< scaledTrace_Type, ExpressionDefinitions::minusTransposedTensor_Type>  scaledTraceTimesMinusTF_Type;
typedef ExpressionProduct< ExpressionScalar,
                           ExpressionProduct<
                               ExpressionDefinitions::deformationGradient_Type,
                               ExpressionDefinitions::outerProduct_Type> >                       scaledFtimesM_Type;
typedef ExpressionProduct< scaledFourthInvariant_Type,
                           ExpressionDefinitions::minusTransposedTensor_Type>                    scaledFourthInvariantTimesMinusTF_Type;

typedef ExpressionAddition<
    ExpressionAddition< scaledTensorF_Type, scaledTraceTimesMinusTF_Type >,
    ExpressionAddition< scaledFtimesM_Type, scaledFourthInvariantTimesMinusTF_Type> > tensorialPart_distrType;
//====================================================
//@}

// Constructors for the expressions defined by the typedefs
//@{

  distributedIsochoricTrace_Type distributedIsochoricTrace( const Real coeff, const ExpressionDefinitions::isochoricTrace_Type ICbar );

  distributedIsochoricStretch_Type distributedIsochoricFourthInvariant( const Real coeff, const ExpressionDefinitions::isochoricStretch_Type I4bar);

 distributedInvariants_Type distributeInvariants( const distributedIsochoricTrace_Type distrIC,
						  const distributedIsochoricStretch_Type distrI4 );

 distributedStretch_Type distributedStretch( const ExpressionDefinitions::isochoricTrace_Type trCBar,
					     const ExpressionDefinitions::isochoricStretch_Type I_4ith, const Real kappa );

  minusTFscalarDF_distrType minusTFscalarDF( const ExpressionDefinitions::minusTransposedTensor_Type minusFT);

  FscalarDF_distrType FscalarDF( const ExpressionDefinitions::deformationGradient_Type F );

 dFTtimesFscalarM_distrType dFTtimesFscalarM( const ExpressionDefinitions::deformationGradient_Type F,
					      const ExpressionDefinitions::outerProduct_Type M);

 FTtimesDFscalarM_distrType FTtimesDFscalarM( const ExpressionDefinitions::deformationGradient_Type F,
					      const ExpressionDefinitions::outerProduct_Type M);

  scaledTrace_Type scaleTrace( const Real coeff, const  ExpressionDefinitions::traceTensor_Type tr );

  scaledIsochoricTrace_Type scaleIsochoricTrace( const Real coeff, const  ExpressionDefinitions::isochoricTrace_Type isoTr );

  scaledDeterminant_Type scaleDeterminant( const Real coeff, const ExpressionDefinitions::powerExpression_Type Jel );

  scaledFourthInvariant_Type scaleFourthInvariant( const Real coeff, const ExpressionDefinitions::stretch_Type I4 );

  scaledIsochoricFourthInvariant_Type scaleIsochoricFourthInvariant( const Real coeff, const ExpressionDefinitions::isochoricStretch_Type isoI4 );

 linearizationFisochoricTrace_Type derivativeIsochoricTrace( const ExpressionDefinitions::isochoricTrace_Type isoTr,
							     const ExpressionDefinitions::powerExpression_Type Jel,
							     const ExpressionDefinitions::deformationGradient_Type F,
							     const ExpressionDefinitions::minusTransposedTensor_Type F_T);

 linearizationFisochoricFourthInvariant_Type derivativeIsochoricFourthInvariant( const ExpressionDefinitions::isochoricStretch_Type isoI4,
										 const ExpressionDefinitions::powerExpression_Type Jel,
										 const ExpressionDefinitions::deformationGradient_Type F,
										 const ExpressionDefinitions::minusTransposedTensor_Type F_T,
										 const ExpressionDefinitions::outerProduct_Type M);

  linearizationDistributedStretch_Type derivativeDistributedStretch( const Real kappa,
								     const ExpressionDefinitions::isochoricTrace_Type isoTr,
								     const ExpressionDefinitions::isochoricStretch_Type isoI4,
								     const ExpressionDefinitions::powerExpression_Type Jel,
								     const ExpressionDefinitions::deformationGradient_Type F,
								     const ExpressionDefinitions::minusTransposedTensor_Type F_T,
								     const ExpressionDefinitions::outerProduct_Type M);

  scaledTensorF_Type scaleF( const Real coeff, const ExpressionDefinitions::deformationGradient_Type F );

 scaledTraceTimesMinusTF_Type scaleTraceMinuTF( const Real coeff,
						const ExpressionDefinitions::traceTensor_Type Ic,
						const ExpressionDefinitions::minusTransposedTensor_Type F_T );

 scaledFtimesM_Type scaleFtimesM( const Real coeff,
                                 const ExpressionDefinitions::deformationGradient_Type F,
				  const ExpressionDefinitions::outerProduct_Type M);

 scaledFourthInvariantTimesMinusTF_Type scaleI4timesMinutTF( const Real coeff,
                                                            const ExpressionDefinitions::stretch_Type I4,
							     const ExpressionDefinitions::minusTransposedTensor_Type F_T );

 tensorialPart_distrType tensorialPartPiola( const Real kappa,
                                            const ExpressionDefinitions::traceTensor_Type tr,
                                            const ExpressionDefinitions::stretch_Type I4,
                                            const ExpressionDefinitions::deformationGradient_Type F,
                                            const ExpressionDefinitions::minusTransposedTensor_Type F_T,
					     const ExpressionDefinitions::outerProduct_Type M);
//@}

} // end namespace ExpressionDistributedModel

// name space that is specific for the Anisotropic Multi-mechanism class.
namespace ExpressionMultimechanism
{
using namespace ExpressionAssembly;

typedef ExpressionSubstraction<
  ExpressionDefinitions::isochoricStretch_Type,
  ExpressionScalar> difference_Type;

typedef ExpressionSubstraction<
    ExpressionDefinitions::stretch_Type,
    ExpressionScalar> incompressibleDifference_Type;

typedef ExpressionDivision<
  incompressibleDifference_Type,
  ExpressionScalar> relativeDifference_Type;

typedef  ExpressionArcTan<difference_Type> activation_Type;

typedef ExpressionVectorFromNonConstantScalar< difference_Type, 3>  expressionVectorFromDifference_Type;
typedef ExpressionVectorFromNonConstantScalar< incompressibleDifference_Type, 3>  expressionVectorFromIncompressibleDifference_Type;
typedef ExpressionVectorFromNonConstantScalar< relativeDifference_Type, 3>  expressionVectorFromRelativeDifference_Type;

typedef ExpressionProduct< ExpressionDefinitions::deformationGradient_Type,
                           ExpressionDefinitions::inverseTensor_Type >   deformationActivatedTensor_Type;

typedef ExpressionProduct< ExpressionDefinitions::minusTransposedTensor_Type,
                           ExpressionTranspose<ExpressionDefinitions::deformationGradient_Type> >   activeMinusTtensor_Type;


typedef ExpressionProduct<
    ExpressionDefinitions::minusTransposedTensor_Type,
    ExpressionProduct< ExpressionDefinitions::rightCauchyGreenTensor_Type,
                       ExpressionDefinitions::inverseTensor_Type > > rightCauchyGreenMultiMechanism_Type;

typedef ExpressionProduct< ExpressionDefinitions::deformationGradient_Type,
                           ExpressionDefinitions::interpolatedValue_Type> activatedFiber_Type;

typedef ExpressionDot< activatedFiber_Type, activatedFiber_Type > squaredNormActivatedFiber_Type;

typedef ExpressionSquareRoot< squaredNormActivatedFiber_Type> normActivatedFiber_Type;

typedef ExpressionNormalize<activatedFiber_Type> normalizedVector_Type;

typedef ExpressionDivision< activatedFiber_Type, normActivatedFiber_Type> normalizedFiber_Type;

typedef ExpressionDivision< ExpressionDefinitions::determinantTensorF_Type,
			    ExpressionDefinitions::interpolatedScalarValue_Type> activatedDeterminantF_Type;

typedef ExpressionProduct< ExpressionDefinitions::determinantTensorF_Type,
			   ExpressionDefinitions::powerExpression_Type> activatedJ_Type;

typedef ExpressionPower<activatedDeterminantF_Type >  activePowerExpression_Type;

typedef ExpressionIsochoricChangeOfVariable<activatedDeterminantF_Type >  activeIsochoricDeterminant_Type;

typedef ExpressionOuterProduct< activatedFiber_Type, activatedFiber_Type>  activeOuterProduct_Type;

typedef ExpressionOuterProduct< normalizedVector_Type, normalizedVector_Type>  activeNormalizedOuterProduct_Type;

typedef ExpressionDot< rightCauchyGreenMultiMechanism_Type, activeNormalizedOuterProduct_Type>  activeStretch_Type;

  typedef ExpressionDot< rightCauchyGreenMultiMechanism_Type, ExpressionDefinitions::outerProduct_Type>  activeInterpolatedFiberStretch_Type;

typedef ExpressionProduct< activeIsochoricDeterminant_Type,
			   activeInterpolatedFiberStretch_Type>         activeIsochoricStretch_Type;

typedef ExpressionProduct< activeIsochoricDeterminant_Type,
			   activeStretch_Type>                          activeNoInterpolationStretch_Type;

typedef ExpressionProduct< activePowerExpression_Type, activeStretch_Type>         activePowerIsochoricStretch_Type;

typedef ExpressionProduct< ExpressionDefinitions::deformationGradient_Type,
                           ExpressionDefinitions::interpolatedValue_Type> activatedFiber_Type;

typedef ExpressionProduct< ExpressionDphiJ,
                           ExpressionDefinitions::inverseTensor_Type>     activeLinearization_Type;


typedef ExpressionProduct< ExpressionDphiI,
                           ExpressionDefinitions::inverseTensor_Type>    activeTestGradient_Type;


difference_Type absoluteStretch( const ExpressionDefinitions::isochoricStretch_Type IVbar,
                                 const Real valueToSubtract);

incompressibleDifference_Type incompressibleAbsoluteStretch(const ExpressionDefinitions::stretch_Type IV,
                                                            const Real valueToSubtract );

relativeDifference_Type relativeDifference(const  incompressibleDifference_Type difference,
					   const Real refFourthInvariant );


activation_Type activationConstructor( const ExpressionMultimechanism::difference_Type absoluteStretch,
                                       const Real intCoeff,
                                       const Real extCoeff,
                                       const Real translation);

expressionVectorFromDifference_Type vectorFromActivation( const ExpressionMultimechanism::difference_Type activation );

expressionVectorFromIncompressibleDifference_Type vectorFromIncompressibleDifference(const
                                                                                     ExpressionMultimechanism::incompressibleDifference_Type activation );

expressionVectorFromRelativeDifference_Type vectorFromRelativeDifference(const
									 ExpressionMultimechanism::relativeDifference_Type relDifference );

deformationActivatedTensor_Type createDeformationActivationTensor( const ExpressionDefinitions::deformationGradient_Type Ft,
                                                                   const ExpressionDefinitions::inverseTensor_Type       F0_ta);

rightCauchyGreenMultiMechanism_Type activationRightCauchyGreen( const ExpressionDefinitions::minusTransposedTensor_Type FzeroAminusT,
                                                                const ExpressionDefinitions::rightCauchyGreenTensor_Type C,
                                                                const ExpressionDefinitions::inverseTensor_Type FzeroAminus1 );

 activatedFiber_Type activateFiberDirection( const ExpressionDefinitions::deformationGradient_Type F,
					     const ExpressionDefinitions::interpolatedValue_Type ithFiber);

 squaredNormActivatedFiber_Type squaredNormActivatedFiber( const activatedFiber_Type f);

 normalizedVector_Type unitVector( const activatedFiber_Type vector );

 normActivatedFiber_Type normActivatedFiber( const activatedFiber_Type f);

 normalizedFiber_Type normalizedFiberDirection( const activatedFiber_Type fiber,
						 const normActivatedFiber_Type normFiber);

 activatedDeterminantF_Type activateDeterminantF( const ExpressionDefinitions::determinantTensorF_Type Jzero,
						  const ExpressionDefinitions::interpolatedScalarValue_Type  JzeroA );

 activatedJ_Type activateJ( const ExpressionDefinitions::determinantTensorF_Type Jzero,
			    const ExpressionDefinitions::powerExpression_Type  JzeroA );

 activePowerExpression_Type activePowerExpression( activatedDeterminantF_Type Ja,
						   const Real exp);

 activeIsochoricDeterminant_Type activeIsochoricDeterminant( activatedDeterminantF_Type Ja );


 activeOuterProduct_Type activeOuterProduct( const activatedFiber_Type activatedFiber );

 activeNormalizedOuterProduct_Type activeNormalizedOuterProduct( const normalizedVector_Type activatedFiber );

 activeStretch_Type activeFiberStretch( const rightCauchyGreenMultiMechanism_Type activeC,
                                        const activeNormalizedOuterProduct_Type activeM);

 activeInterpolatedFiberStretch_Type activeInterpolatedFiberStretch( const rightCauchyGreenMultiMechanism_Type activeC,
								const ExpressionDefinitions::outerProduct_Type activeM);

 activeIsochoricStretch_Type activeIsochoricFourthInvariant( const activeIsochoricDeterminant_Type activeJ,
							     const activeInterpolatedFiberStretch_Type activeI4);

 activeNoInterpolationStretch_Type activeNoInterpolationFourthInvariant( const activeIsochoricDeterminant_Type activeJ,
									 const activeStretch_Type activeI4);

 activePowerIsochoricStretch_Type activePowerIsochoricFourthInvariant( const activePowerExpression_Type activeJ,
								       const activeStretch_Type activeI4);

 activeMinusTtensor_Type createActiveMinusTtensor( const ExpressionDefinitions::minusTransposedTensor_Type FminusT,
						   const ExpressionTranspose<ExpressionDefinitions::deformationGradient_Type> FzeroA);

 activeLinearization_Type activatedLinearization( const ExpressionDphiJ der,
						  const ExpressionDefinitions::inverseTensor_Type inverse);

 activeTestGradient_Type activatedTestGradient(const ExpressionDphiI gradTest,
					       const ExpressionDefinitions::inverseTensor_Type FAminus1);
}// end namespace ExpressionDistributedModel

} //! End namespace LifeV
#endif
