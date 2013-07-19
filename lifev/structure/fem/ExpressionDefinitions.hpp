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

// Definition of the isochoric trace \bar{I_C} = J^( -2.0/3.0)*tr(C))
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

namespace ExpressionMultimechanism
{
using namespace ExpressionAssembly;

typedef ExpressionSubstraction<
    ExpressionDefinitions::isochoricStretch_Type,
    ExpressionScalar> difference_Type;

typedef  ExpressionArcTan<difference_Type> activation_Type;

typedef ExpressionVectorFromNonConstantScalar< difference_Type, 3>  expressionVectorFromDifference_Type;

typedef ExpressionProduct< ExpressionDefinitions::deformationGradient_Type,
                           ExpressionDefinitions::inverseTensor_Type >   deformationActivatedTensor_Type;


typedef ExpressionProduct<
    ExpressionDefinitions::minusTransposedTensor_Type,
    ExpressionProduct< ExpressionDefinitions::rightCauchyGreenTensor_Type,
                       ExpressionDefinitions::inverseTensor_Type > > rightCauchyGreenMultiMechanism_Type;

typedef ExpressionProduct< ExpressionDefinitions::deformationGradient_Type,
                           ExpressionDefinitions::interpolatedValue_Type> activatedFiber_Type;


typedef ExpressionProduct< ExpressionDefinitions::determinantTensorF_Type,
                           ExpressionDefinitions::determinantTensorF_Type> activatedDeterminantF_Type;

typedef ExpressionPower<activatedDeterminantF_Type >  activePowerExpression_Type;

typedef ExpressionOuterProduct< activatedFiber_Type, activatedFiber_Type>  activeOuterProduct_Type;

typedef ExpressionDot< rightCauchyGreenMultiMechanism_Type, activeOuterProduct_Type>  activeStretch_Type;

typedef ExpressionProduct< activePowerExpression_Type, activeStretch_Type>         activeIsochoricStretch_Type;

difference_Type absoluteStretch( const ExpressionDefinitions::isochoricStretch_Type IVbar,
                                 const Real valueToSubtract)
{
    return difference_Type ( IVbar, value( valueToSubtract ) );

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

activatedDeterminantF_Type activateDeterminantF( const ExpressionDefinitions::determinantTensorF_Type Jzero,
                                                 const ExpressionDefinitions::determinantTensorF_Type JzeroA )
{
    return activatedDeterminantF_Type( Jzero, JzeroA );
}

activePowerExpression_Type activePowerExpression( activatedDeterminantF_Type Ja,
                                                  const Real exp)
{
    return activePowerExpression_Type ( Ja, exp );
}

activeOuterProduct_Type activeOuterProduct( const activatedFiber_Type activatedFiber )
{
    return activeOuterProduct_Type( activatedFiber, activatedFiber );
}

activeStretch_Type activeFiberStretch( const rightCauchyGreenMultiMechanism_Type activeC,
                                         const activeOuterProduct_Type activeM)
{
    return activeStretch_Type( activeC, activeM );
}

activeIsochoricStretch_Type activeIsochoricFourthInvariant( const activePowerExpression_Type activeJ,
                                                            const activeStretch_Type activeI4)
{
    return activeIsochoricStretch_Type( activeJ, activeI4 );
}

}// end namespace ExpressionDistributedModel
#endif

} //! End namespace LifeV
#endif
