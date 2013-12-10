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


#ifndef _EXPRESSIONDEFINITIONS_CPP_
#define _EXPRESSIONDEFINITIONS_CPP_

#include <lifev/structure/fem/ExpressionDefinitions.hpp>
#include <boost/shared_ptr.hpp>


namespace LifeV
{
  namespace ExpressionDefinitions
  {

    using namespace ExpressionAssembly;

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

    inverseTensor_Type inv( const deformationGradient_Type F )
    {
      return inverseTensor_Type( F );
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

    isochoricChangeOfVariable_Type isochoricDeterminant( const determinantTensorF_Type J )
    {
      return isochoricChangeOfVariable_Type( J );
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

    interpolatedValue_Type interpolateValue( const boost::shared_ptr< ETFESpace_Type > dispETFESpace,
					     const vector_Type& valueVector)
    {
      return interpolatedValue_Type( dispETFESpace, valueVector ) ;
    }

    interpolatedScalarValue_Type interpolateScalarValue( const boost::shared_ptr< scalarETFESpace_Type > dispETFESpace,
							 const vector_Type& valueVector)
    {
      return interpolatedScalarValue_Type( dispETFESpace, valueVector ) ;
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

#endif
  } //! End namespace ExpressionDefinitions


  //! The namespace ExpressionDistributedModel is specific for the Dstributed Holzapfel model
  //! the definitions have been inserted here in order to avoid huge declarations of expressions
  //! in the header file of the model

  //! The goal of the current namespace is just to define the expressions (final and intermediate)
  //! that are needed for the distributed model making use of the previous namespaces already
  //! defined.

#ifdef ENABLE_ANISOTROPIC_LAW
  namespace ExpressionDistributedModel
  {
    using namespace ExpressionAssembly;

    // Constructors for the expressions defined by the typedefs
    //@{

    distributedIsochoricTrace_Type distributedIsochoricTrace( const Real coeff, const ExpressionDefinitions::isochoricTrace_Type ICbar )
    {
      return distributedIsochoricTrace_Type( value( coeff ), ICbar);
    }

    distributedIsochoricStretch_Type distributedIsochoricFourthInvariant( const Real coeff, const ExpressionDefinitions::isochoricStretch_Type I4bar)
    {
      return distributedIsochoricStretch_Type( value( coeff ), I4bar );
    }

    distributedInvariants_Type distributeInvariants( const distributedIsochoricTrace_Type distrIC,
						     const distributedIsochoricStretch_Type distrI4 )
    {
      return distributedInvariants_Type( distrIC,distrI4 );
    }

    distributedStretch_Type distributedStretch( const ExpressionDefinitions::isochoricTrace_Type trCBar,
						const ExpressionDefinitions::isochoricStretch_Type I_4ith, const Real kappa )
    {
      distributedIsochoricTrace_Type dIC_bar = distributedIsochoricTrace( kappa, trCBar ) ;
      distributedIsochoricStretch_Type dI4_bar = distributedIsochoricFourthInvariant( 1.0 - 3.0 * kappa, I_4ith ) ;

      distributedInvariants_Type dInvariants = distributeInvariants( dIC_bar, dI4_bar );

      return distributedStretch_Type( dInvariants, value( -1.0 ) );
    }

    minusTFscalarDF_distrType minusTFscalarDF( const ExpressionDefinitions::minusTransposedTensor_Type minusFT)
    {
      return minusTFscalarDF_distrType( minusFT, grad( phi_j ) );
    }

    FscalarDF_distrType FscalarDF( const ExpressionDefinitions::deformationGradient_Type F )
    {
      return FscalarDF_distrType( F, grad(phi_j) );
    }

    dFTtimesFscalarM_distrType dFTtimesFscalarM( const ExpressionDefinitions::deformationGradient_Type F,
						 const ExpressionDefinitions::outerProduct_Type M)
    {
      return dFTtimesFscalarM_distrType( transpose( grad(phi_j) ) * F, M );
    }

    FTtimesDFscalarM_distrType FTtimesDFscalarM( const ExpressionDefinitions::deformationGradient_Type F,
						 const ExpressionDefinitions::outerProduct_Type M)
    {
      return FTtimesDFscalarM_distrType( transpose( F ) * grad( phi_j ), M );
    }

    scaledTrace_Type scaleTrace( const Real coeff, const  ExpressionDefinitions::traceTensor_Type tr )
    {
      return scaledTrace_Type( value( coeff ), tr );
    }

    scaledIsochoricTrace_Type scaleIsochoricTrace( const Real coeff, const  ExpressionDefinitions::isochoricTrace_Type isoTr )
    {
      return scaledIsochoricTrace_Type( value( coeff ), isoTr );
    }

    scaledDeterminant_Type scaleDeterminant( const Real coeff, const ExpressionDefinitions::powerExpression_Type Jel )
    {
      return scaledDeterminant_Type( value( coeff ), Jel );
    }

    scaledFourthInvariant_Type scaleFourthInvariant( const Real coeff, const ExpressionDefinitions::stretch_Type I4 )
    {
      return scaledFourthInvariant_Type( value( coeff ), I4 );
    }

    scaledIsochoricFourthInvariant_Type scaleIsochoricFourthInvariant( const Real coeff, const ExpressionDefinitions::isochoricStretch_Type isoI4 )
    {
      return scaledIsochoricFourthInvariant_Type( value( coeff ) , isoI4 );
    }

    linearizationFisochoricTrace_Type derivativeIsochoricTrace( const ExpressionDefinitions::isochoricTrace_Type isoTr,
								const ExpressionDefinitions::powerExpression_Type Jel,
								const ExpressionDefinitions::deformationGradient_Type F,
								const ExpressionDefinitions::minusTransposedTensor_Type F_T)
    {
      scaledIsochoricTrace_Type scIsoTr = scaleIsochoricTrace(  ( -2.0/3.0 ), isoTr );
      scaledDeterminant_Type    scJ = scaleDeterminant( 2.0, Jel );

      minusTFscalarDF_distrType FTdotDF = minusTFscalarDF( F_T );
      FscalarDF_distrType       FdotDF = FscalarDF( F );

      return linearizationFisochoricTrace_Type( scIsoTr * FTdotDF , scJ * FdotDF );
    }

    linearizationFisochoricFourthInvariant_Type derivativeIsochoricFourthInvariant( const ExpressionDefinitions::isochoricStretch_Type isoI4,
										    const ExpressionDefinitions::powerExpression_Type Jel,
										    const ExpressionDefinitions::deformationGradient_Type F,
										    const ExpressionDefinitions::minusTransposedTensor_Type F_T,
										    const ExpressionDefinitions::outerProduct_Type M)
    {
      scaledIsochoricFourthInvariant_Type scIsoI4 = scaleIsochoricFourthInvariant( ( -2.0/3.0 ), isoI4 );

      minusTFscalarDF_distrType FTdotDF = minusTFscalarDF( F_T );

      dFTtimesFscalarM_distrType firstTerm = dFTtimesFscalarM( F, M );
      FTtimesDFscalarM_distrType secondTerm = FTtimesDFscalarM( F, M );

      return linearizationFisochoricFourthInvariant_Type( scIsoI4 * FTdotDF , Jel * ( firstTerm + secondTerm ) );
    }

    linearizationDistributedStretch_Type derivativeDistributedStretch( const Real kappa,
								       const ExpressionDefinitions::isochoricTrace_Type isoTr,
								       const ExpressionDefinitions::isochoricStretch_Type isoI4,
								       const ExpressionDefinitions::powerExpression_Type Jel,
								       const ExpressionDefinitions::deformationGradient_Type F,
								       const ExpressionDefinitions::minusTransposedTensor_Type F_T,
								       const ExpressionDefinitions::outerProduct_Type M)
    {
      linearizationFisochoricTrace_Type firstTerm = derivativeIsochoricTrace( isoTr, Jel, F, F_T );
      linearizationFisochoricFourthInvariant_Type secondTerm = derivativeIsochoricFourthInvariant( isoI4, Jel, F, F_T, M);

      return linearizationDistributedStretch_Type( value(kappa) * firstTerm, value( 1.0 - 3.0 * kappa) * secondTerm );
    }

    scaledTensorF_Type scaleF( const Real coeff, const ExpressionDefinitions::deformationGradient_Type F )
    {
      return scaledTensorF_Type( coeff, F);
    }

    scaledTraceTimesMinusTF_Type scaleTraceMinuTF( const Real coeff,
						   const ExpressionDefinitions::traceTensor_Type Ic,
						   const ExpressionDefinitions::minusTransposedTensor_Type F_T )
    {
      scaledTrace_Type scIc = scaleTrace( coeff, Ic );

      return scaledTraceTimesMinusTF_Type( scIc, F_T );
    }

    scaledFtimesM_Type scaleFtimesM( const Real coeff,
				     const ExpressionDefinitions::deformationGradient_Type F,
				     const ExpressionDefinitions::outerProduct_Type M)
    {

      return scaledFtimesM_Type( coeff, F*M );
    }

    scaledFourthInvariantTimesMinusTF_Type scaleI4timesMinutTF( const Real coeff,
								const ExpressionDefinitions::stretch_Type I4,
								const ExpressionDefinitions::minusTransposedTensor_Type F_T )
    {
      scaledFourthInvariant_Type scI4 = scaleFourthInvariant( coeff , I4 );

      return scaledFourthInvariantTimesMinusTF_Type( scI4, F_T );
    }

    tensorialPart_distrType tensorialPartPiola( const Real kappa,
						const ExpressionDefinitions::traceTensor_Type tr,
						const ExpressionDefinitions::stretch_Type I4,
						const ExpressionDefinitions::deformationGradient_Type F,
						const ExpressionDefinitions::minusTransposedTensor_Type F_T,
						const ExpressionDefinitions::outerProduct_Type M)
    {
      /*
	First the terms of the expression are built term by term
	The tensorial part of the Piola tensor corresponding to one fiber
	reads:

	kappa * F - ( kappa / 3.0 ) * I_C * F_T + ( 1 - 3*kappa) FM - ( 1 - 3*kappa )/3.0 * I4 * F_T
      */

      scaledTensorF_Type scF = scaleF( kappa, F );
      scaledTraceTimesMinusTF_Type scTrF_T = scaleTraceMinuTF( - kappa / 3.0, tr, F_T );

      scaledFtimesM_Type scFM = scaleFtimesM( ( 1 - 3.0 * kappa ), F, M );
      scaledFourthInvariantTimesMinusTF_Type scI4F_T = scaleI4timesMinutTF( - ( 1.0 - 3.0 * kappa ) / 3.0, I4, F_T );

      return tensorialPart_distrType( ( scF + scTrF_T ), ( scFM + scI4F_T ) );
    }
    //@}

  } // end namespace ExpressionDistributedModel

  // name space that is specific for the Anisotropic Multi-mechanism class.
  namespace ExpressionMultimechanism
  {
    using namespace ExpressionAssembly;

    difference_Type absoluteStretch( const ExpressionDefinitions::isochoricStretch_Type IVbar,
				     const Real valueToSubtract)
    {
      return difference_Type ( IVbar, value( valueToSubtract ) );

    }

    incompressibleDifference_Type incompressibleAbsoluteStretch( const ExpressionDefinitions::stretch_Type IV,
                                                                 const Real valueToSubtract)
    {
      return incompressibleDifference_Type ( IV, value( valueToSubtract ) );

    }

    relativeDifference_Type relativeDifference( const  incompressibleDifference_Type difference,
						const Real refFourthInvariant)
    {
      return relativeDifference_Type ( difference, value( refFourthInvariant ) );

    }


    activation_Type activationConstructor( const ExpressionMultimechanism::difference_Type absoluteStretch,
					   const Real intCoeff,
					   const Real extCoeff,
					   const Real translation)
    {
      return activation_Type( absoluteStretch, intCoeff, extCoeff, translation );
    }

    expressionVectorFromDifference_Type vectorFromActivation( const ExpressionMultimechanism::difference_Type activation )
    {
      return expressionVectorFromDifference_Type( activation );
    }

    expressionVectorFromIncompressibleDifference_Type vectorFromIncompressibleDifference( const ExpressionMultimechanism::incompressibleDifference_Type activation )
    {
      return expressionVectorFromIncompressibleDifference_Type( activation );
    }

    expressionVectorFromRelativeDifference_Type vectorFromRelativeDifference( const ExpressionMultimechanism::relativeDifference_Type activation )
    {
      return expressionVectorFromRelativeDifference_Type( activation );
    }

    deformationActivatedTensor_Type createDeformationActivationTensor( const ExpressionDefinitions::deformationGradient_Type Ft,
								       const ExpressionDefinitions::inverseTensor_Type       F0_ta)
    {
      return deformationActivatedTensor_Type( Ft, F0_ta );
    }

    rightCauchyGreenMultiMechanism_Type activationRightCauchyGreen( const ExpressionDefinitions::minusTransposedTensor_Type FzeroAminusT,
								    const ExpressionDefinitions::rightCauchyGreenTensor_Type C,
								    const ExpressionDefinitions::inverseTensor_Type FzeroAminus1 )
    {
      ExpressionProduct< ExpressionDefinitions::rightCauchyGreenTensor_Type,
			 ExpressionDefinitions::inverseTensor_Type > rightTerm( C, FzeroAminus1 );
      return rightCauchyGreenMultiMechanism_Type( FzeroAminusT, rightTerm );
    }

    activatedFiber_Type activateFiberDirection( const ExpressionDefinitions::deformationGradient_Type F,
						const ExpressionDefinitions::interpolatedValue_Type ithFiber)
    {
      return activatedFiber_Type( F, ithFiber );
    }

    normalizedVector_Type unitVector( const activatedFiber_Type vector)
    {
      return normalizedVector_Type( vector );
    }

    squaredNormActivatedFiber_Type squaredNormActivatedFiber( const activatedFiber_Type f)
    {
      return squaredNormActivatedFiber_Type( f, f );
    }

    normActivatedFiber_Type normActivatedFiber( const activatedFiber_Type f)
    {
      squaredNormActivatedFiber_Type squaredNorm = squaredNormActivatedFiber_Type( f, f );

      return normActivatedFiber_Type( squaredNorm );
    }

    normalizedFiber_Type normalizedFiberDirection( const activatedFiber_Type fiber,
						   const normActivatedFiber_Type normFiber)
    {
      return normalizedFiber_Type( fiber, normFiber );
    }

    activatedDeterminantF_Type activateDeterminantF( const ExpressionDefinitions::determinantTensorF_Type Jzero,
						     const ExpressionDefinitions::interpolatedScalarValue_Type  JzeroA )
    {
      return activatedDeterminantF_Type( Jzero, JzeroA );
    }

    activatedJ_Type activateJ( const ExpressionDefinitions::determinantTensorF_Type Jzero,
			       const ExpressionDefinitions::powerExpression_Type  JzeroA )
    {
      return activatedJ_Type( Jzero, JzeroA );
    }

    activePowerExpression_Type activePowerExpression( activatedDeterminantF_Type Ja,
						      const Real exp)
    {
      return activePowerExpression_Type ( Ja, exp );
    }

    activeIsochoricDeterminant_Type activeIsochoricDeterminant( activatedDeterminantF_Type Ja )
    {
      return activeIsochoricDeterminant_Type ( Ja );
    }

    activeOuterProduct_Type activeOuterProduct( const activatedFiber_Type activatedFiber )
    {
      return activeOuterProduct_Type( activatedFiber, activatedFiber );
    }

    activeNormalizedOuterProduct_Type activeNormalizedOuterProduct( const normalizedVector_Type normalizedActiveFiber )
    {
      return activeNormalizedOuterProduct_Type( normalizedActiveFiber, normalizedActiveFiber );
    }

    activeStretch_Type activeFiberStretch( const rightCauchyGreenMultiMechanism_Type activeC,
                                           const activeNormalizedOuterProduct_Type activeM)
    {
      return activeStretch_Type( activeC, activeM );
    }

    activeInterpolatedFiberStretch_Type activeInterpolatedFiberStretch( const rightCauchyGreenMultiMechanism_Type activeC,
									const ExpressionDefinitions::outerProduct_Type activeM)
    {
      return activeInterpolatedFiberStretch_Type( activeC, activeM );
    }

    activeIsochoricStretch_Type activeIsochoricFourthInvariant( const activeIsochoricDeterminant_Type activeJ,
								const activeInterpolatedFiberStretch_Type activeI4)
    {
      return activeIsochoricStretch_Type( activeJ, activeI4 );
    }

    activeNoInterpolationStretch_Type activeNoInterpolationFourthInvariant( const activeIsochoricDeterminant_Type activeJ,
									    const activeStretch_Type activeI4)
    {
      return activeNoInterpolationStretch_Type( activeJ, activeI4 );
    }

    activePowerIsochoricStretch_Type activePowerIsochoricFourthInvariant( const activePowerExpression_Type activeJ,
									   const activeStretch_Type activeI4)
    {
      return activePowerIsochoricStretch_Type( activeJ, activeI4 );
    }


    activeMinusTtensor_Type createActiveMinusTtensor( const ExpressionDefinitions::minusTransposedTensor_Type FminusT,
						      const ExpressionTranspose<ExpressionDefinitions::deformationGradient_Type> FzeroA)
    {
      return activeMinusTtensor_Type ( FminusT, FzeroA );
    }

    activeLinearization_Type activatedLinearization( const ExpressionDphiJ der,
						     const ExpressionDefinitions::inverseTensor_Type inverse)
    {
      return activeLinearization_Type( der , inverse);
    }
    activeTestGradient_Type activatedTestGradient(const ExpressionDphiI gradTest,
						  const ExpressionDefinitions::inverseTensor_Type FAminus1)
    {
      return activeTestGradient_Type ( gradTest, FAminus1 );
    }
  }// end namespace ExpressionDistributedModel
#endif

} //! End namespace LifeV
#endif
