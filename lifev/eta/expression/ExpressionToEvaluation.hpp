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
     @brief This file contains the definition of the ExpressionToEvaluation class

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef EXPRESSION_TO_EVALUATION_HPP
#define EXPRESSION_TO_EVALUATION_HPP

#include <lifev/eta/expression/ExpressionPhiI.hpp>
#include <lifev/eta/expression/ExpressionPhiJ.hpp>
#include <lifev/eta/expression/ExpressionDphiI.hpp>
#include <lifev/eta/expression/ExpressionDphiJ.hpp>
#include <lifev/eta/expression/ExpressionDivI.hpp>
#include <lifev/eta/expression/ExpressionDivJ.hpp>
#include <lifev/eta/expression/ExpressionMinusTransposed.hpp>
#include <lifev/eta/expression/ExpressionInverse.hpp>
#include <lifev/eta/expression/ExpressionDeterminant.hpp>
#include <lifev/eta/expression/ExpressionTrace.hpp>
#include <lifev/eta/expression/ExpressionLaplacianI.hpp>
#include <lifev/eta/expression/ExpressionLaplacianJ.hpp>

#include <lifev/eta/expression/ExpressionAddition.hpp>
#include <lifev/eta/expression/ExpressionSubstraction.hpp>
#include <lifev/eta/expression/ExpressionProduct.hpp>
#include <lifev/eta/expression/ExpressionPower.hpp>
#include <lifev/eta/expression/ExpressionSquareRoot.hpp>
#include <lifev/eta/expression/ExpressionCubicRoot.hpp>
#include <lifev/eta/expression/ExpressionArcTan.hpp>
#include <lifev/eta/expression/ExpressionDerivativeArcTan.hpp>
#include <lifev/eta/expression/ExpressionLogarithm.hpp>
#include <lifev/eta/expression/ExpressionExponential.hpp>
#include <lifev/eta/expression/ExpressionDot.hpp>
#include <lifev/eta/expression/ExpressionDivision.hpp>
#include <lifev/eta/expression/ExpressionEmult.hpp>
#include <lifev/eta/expression/ExpressionExtract1.hpp>
#include <lifev/eta/expression/ExpressionExtract2.hpp>
#include <lifev/eta/expression/ExpressionTranspose.hpp>
#include <lifev/eta/expression/ExpressionNormalize.hpp>
#include <lifev/eta/expression/ExpressionSymmetricTensor.hpp>
#include <lifev/eta/expression/ExpressionOuterProduct.hpp>

#include <lifev/eta/expression/ExpressionScalar.hpp>
#include <lifev/eta/expression/ExpressionVector.hpp>
#include <lifev/eta/expression/ExpressionMatrix.hpp>

#include <lifev/eta/expression/ExpressionInterpolateValue.hpp>
#include <lifev/eta/expression/ExpressionInterpolateGradient.hpp>
#include <lifev/eta/expression/ExpressionInterpolateLaplacian.hpp>

#include <lifev/eta/expression/ExpressionFunctor.hpp>

#include <lifev/eta/expression/ExpressionHK.hpp>
#include <lifev/eta/expression/ExpressionMeas.hpp>
#include <lifev/eta/expression/ExpressionMeasBDCurrentFE.hpp>
#include <lifev/eta/expression/ExpressionPosition.hpp>
#include <lifev/eta/expression/ExpressionNormal.hpp>

#include <lifev/eta/expression/ExpressionIfCrossed.hpp>
#include <lifev/eta/expression/ExpressionVectorFromNonConstantScalar.hpp>
#include <lifev/eta/expression/ExpressionVectorFromNonConstantMatrix.hpp>
#include <lifev/eta/expression/ExpressionPatchArea.hpp>

#include <lifev/eta/expression/ExpressionScalarToVector.hpp>

#include <lifev/eta/expression/EvaluationPhiI.hpp>
#include <lifev/eta/expression/EvaluationPhiJ.hpp>
#include <lifev/eta/expression/EvaluationDphiI.hpp>
#include <lifev/eta/expression/EvaluationDphiJ.hpp>
#include <lifev/eta/expression/EvaluationDivI.hpp>
#include <lifev/eta/expression/EvaluationDivJ.hpp>
#include <lifev/eta/expression/EvaluationMinusTransposed.hpp>
#include <lifev/eta/expression/EvaluationInverse.hpp>
#include <lifev/eta/expression/EvaluationDeterminant.hpp>
#include <lifev/eta/expression/EvaluationTrace.hpp>
#include <lifev/eta/expression/EvaluationLaplacianI.hpp>
#include <lifev/eta/expression/EvaluationLaplacianJ.hpp>

#include <lifev/eta/expression/EvaluationAddition.hpp>
#include <lifev/eta/expression/EvaluationSubstraction.hpp>
#include <lifev/eta/expression/EvaluationProduct.hpp>
#include <lifev/eta/expression/EvaluationPower.hpp>
#include <lifev/eta/expression/EvaluationSquareRoot.hpp>
#include <lifev/eta/expression/EvaluationCubicRoot.hpp>
#include <lifev/eta/expression/EvaluationArcTan.hpp>
#include <lifev/eta/expression/EvaluationDerivativeArcTan.hpp>
#include <lifev/eta/expression/EvaluationLogarithm.hpp>
#include <lifev/eta/expression/EvaluationExponential.hpp>
#include <lifev/eta/expression/EvaluationDot.hpp>
#include <lifev/eta/expression/EvaluationDivision.hpp>
#include <lifev/eta/expression/EvaluationEmult.hpp>
#include <lifev/eta/expression/EvaluationExtract1.hpp>
#include <lifev/eta/expression/EvaluationExtract2.hpp>
#include <lifev/eta/expression/EvaluationTranspose.hpp>
#include <lifev/eta/expression/EvaluationNormalize.hpp>
#include <lifev/eta/expression/EvaluationSymmetricTensor.hpp>
#include <lifev/eta/expression/EvaluationOuterProduct.hpp>

#include <lifev/eta/expression/EvaluationScalar.hpp>
#include <lifev/eta/expression/EvaluationVector.hpp>
#include <lifev/eta/expression/EvaluationMatrix.hpp>

#include <lifev/eta/expression/EvaluationInterpolateValue.hpp>
#include <lifev/eta/expression/EvaluationInterpolateGradient.hpp>
#include <lifev/eta/expression/EvaluationInterpolateLaplacian.hpp>

#include <lifev/eta/expression/EvaluationFunctor.hpp>

#include <lifev/eta/expression/EvaluationHK.hpp>
#include <lifev/eta/expression/EvaluationMeas.hpp>
#include <lifev/eta/expression/EvaluationMeasBDCurrentFE.hpp>
#include <lifev/eta/expression/EvaluationPosition.hpp>
#include <lifev/eta/expression/EvaluationNormal.hpp>

#include <lifev/eta/expression/EvaluationIfCrossed.hpp>
#include <lifev/eta/expression/EvaluationPatchArea.hpp>
#include <lifev/eta/expression/EvaluationVectorFromNonConstantScalar.hpp>
#include <lifev/eta/expression/EvaluationVectorFromNonConstantMatrix.hpp>

#include <lifev/eta/expression/EvaluationScalarToVector.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionToEvaluation  A class to pass from an Expression (Tree) to the corresponding Evaluation (Tree)
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class cannot be instanciated, neither its generic definition nor its partial or total specializations.

  The only scope of this class is to provide, given an Expression (Tree), the type of the corresponding
  Evaluation (Tree). This is achieved by a typedef "evaluation_Type" provided only in the partial specializations.

  The template parameters are the following (no requirements):

  <b> Expression </b> The type of the expression to be "converted" to an Evaluation
  <b> testDim </b> The field dimension of the test finite element space
  <b> solutionDim </b> The field dimension of the solution finite element space
  <b> spaceDim </b> The dimension of the domain.

  Remark that the field dimension of the test finite element space and the field dimension of the solution finite element
  space can differ. This is the case e.g. when assembling a pressure/velocity block for a Stokes problem.

*/

template<typename Expression, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation
{
public:
    typedef typename Expression::evaluation_Type evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// \cond

// Specialized for Extract1
template<typename Expression, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionExtract1<Expression>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationExtract1 <
    typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for Extract2
template<typename Expression, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionExtract2<Expression>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationExtract2 <
    typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for phi_i
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionPhiI, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationPhiI<testDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for phi_j
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionPhiJ, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationPhiJ<solutionDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for dphi_i
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionDphiI, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationDphiI<testDim, spaceDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for laplacian(phi_i)
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionLaplacianI, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationLaplacianI<testDim, spaceDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for laplacian(phi_j)
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionLaplacianJ, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationLaplacianJ<solutionDim, spaceDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for dphi_j
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionDphiJ, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationDphiJ<solutionDim, spaceDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for div_i
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionDivI, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationDivI<testDim, spaceDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for div_j
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionDivJ, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationDivJ<solutionDim, spaceDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for transpose
template<typename Expression, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionTranspose<Expression>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationTranspose <
    typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for normalization
template<typename Expression, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionNormalize<Expression>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationNormalize <
    typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for symmetric expression
template<typename Expression, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionSymmetricTensor<Expression>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationSymmetricTensor <
    typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for the minus transposed
template<typename Expression, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionMinusTransposed<Expression>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationMinusTransposed <
    typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};
    
// 
template<typename Expression, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
ExpressionScalarToVector<Expression, spaceDim>
, testDim
, solutionDim
, spaceDim >
{
public:
    typedef ExpressionScalarToVector <
    typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type
    , spaceDim > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for the inverse
template<typename Expression, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionInverse<Expression>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationInverse <
    typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for determinant
template<typename Expression, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionDeterminant<Expression>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationDeterminant <
    typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for trace
template<typename Expression, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionTrace<Expression>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationTrace <
    typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for vector from non constant scalar fields
template<typename Expression, UInt FEFieldDim, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
  ExpressionVectorFromNonConstantScalar<Expression, FEFieldDim>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationVectorFromNonConstantScalar <
  typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type,
  FEFieldDim
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for vector from non constant scalar fields
template<typename Expression, UInt FESpaceDim, UInt FEFieldDim, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
  ExpressionVectorFromNonConstantMatrix<Expression, FESpaceDim, FEFieldDim>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationVectorFromNonConstantMatrix<
  typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type,
  FESpaceDim,
  FEFieldDim
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for scalar
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionScalar, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationScalar evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for scalar given by a functor
template<typename VectorType, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionExtractScalar<VectorType>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationExtractScalar<VectorType>    evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for vector
template<UInt testDim, UInt solutionDim, UInt spaceDim, UInt VectorDim>
class ExpressionToEvaluation<ExpressionVector<VectorDim>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationVector<VectorDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for matrix
template<UInt testDim, UInt solutionDim, UInt spaceDim, UInt MatrixDim1, UInt MatrixDim2>
class ExpressionToEvaluation<ExpressionMatrix<MatrixDim1, MatrixDim2>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationMatrix<MatrixDim1, MatrixDim2> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for an interpolated value
template<typename MeshType, typename MapType, UInt FESpaceDim, UInt FEFieldDim, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionInterpolateValue<MeshType, MapType, FESpaceDim, FEFieldDim>, testDim, solutionDim, spaceDim >
{
public:
    typedef EvaluationInterpolateValue<MeshType, MapType, FESpaceDim, FEFieldDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for a patch area expression
template<typename MeshType, typename MapType, UInt FESpaceDim, UInt FEFieldDim, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionPatchArea<MeshType, MapType, FESpaceDim, FEFieldDim>, testDim, solutionDim, spaceDim >
{
public:
    typedef EvaluationPatchArea<MeshType, MapType, FESpaceDim, FEFieldDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for an interpolated gradient
template<typename MeshType, typename MapType, UInt FESpaceDim, UInt FEFieldDim, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionInterpolateGradient<MeshType, MapType, FESpaceDim, FEFieldDim>, testDim, solutionDim, spaceDim >
{
public:
    typedef EvaluationInterpolateGradient<MeshType, MapType, FESpaceDim, FEFieldDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for an interpolated laplacian
template<typename MeshType, typename MapType, UInt FESpaceDim, UInt FEFieldDim, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionInterpolateLaplacian<MeshType, MapType, FESpaceDim, FEFieldDim>, testDim, solutionDim, spaceDim >
{
public:
    typedef EvaluationInterpolateLaplacian<MeshType, MapType, FESpaceDim, FEFieldDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for a functor
template<typename Functor, typename Expression, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionFunctor1<Functor, Expression>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationFunctor1 <
    Functor,
    typename ExpressionToEvaluation<Expression, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for a functor 2 arguments
template<typename Functor, typename Expression1, typename Expression2, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionFunctor2<Functor, Expression1, Expression2>
    , testDim
    , solutionDim
    , spaceDim >
{
public:
    typedef EvaluationFunctor2 <
    Functor,
    typename ExpressionToEvaluation<Expression1, testDim, solutionDim, spaceDim>::evaluation_Type,
    typename ExpressionToEvaluation<Expression2, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for a sum
template<typename ExpressionL, typename ExpressionR, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionAddition<ExpressionL, ExpressionR>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationAddition <
    typename ExpressionToEvaluation<ExpressionL, testDim, solutionDim, spaceDim>::evaluation_Type
    , typename ExpressionToEvaluation<ExpressionR, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for a substraction
template<typename ExpressionL, typename ExpressionR, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionSubstraction<ExpressionL, ExpressionR>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationSubstraction <
    typename ExpressionToEvaluation<ExpressionL, testDim, solutionDim, spaceDim>::evaluation_Type
    , typename ExpressionToEvaluation<ExpressionR, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for a product
template<typename ExpressionL, typename ExpressionR, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionProduct<ExpressionL, ExpressionR>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationProduct <
    typename ExpressionToEvaluation<ExpressionL, testDim, solutionDim, spaceDim>::evaluation_Type
    , typename ExpressionToEvaluation<ExpressionR, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for a power
template<typename ExpressionBase, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionPower<ExpressionBase>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationPower <
    typename ExpressionToEvaluation<ExpressionBase, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for the square root
template<typename ExpressionBase, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionSquareRoot<ExpressionBase>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationSquareRoot <
    typename ExpressionToEvaluation<ExpressionBase, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for the cubic root
template<typename ExpressionBase, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionCubicRoot<ExpressionBase>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationCubicRoot <
    typename ExpressionToEvaluation<ExpressionBase, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for the isochoric change of variable in the structure module
template<typename ExpressionBase, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionIsochoricChangeOfVariable<ExpressionBase>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationIsochoricChangeOfVariable <
    typename ExpressionToEvaluation<ExpressionBase, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for the atan
template<typename ExpressionBase, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionArcTan<ExpressionBase>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationArcTan <
    typename ExpressionToEvaluation<ExpressionBase, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for the atan
template<typename ExpressionBase, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionDerivativeArcTan<ExpressionBase>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationDerivativeArcTan <
    typename ExpressionToEvaluation<ExpressionBase, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};


// Specialized for a logarithm
template<typename ExpressionBase, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionLogarithm<ExpressionBase>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationLogarithm <
    typename ExpressionToEvaluation<ExpressionBase, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for a exponential
template<typename ExpressionBase, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionExponential<ExpressionBase>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationExponential <
    typename ExpressionToEvaluation<ExpressionBase, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for a dot product
template<typename ExpressionL, typename ExpressionR, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionDot<ExpressionL, ExpressionR>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationDot <
    typename ExpressionToEvaluation<ExpressionL, testDim, solutionDim, spaceDim>::evaluation_Type
    , typename ExpressionToEvaluation<ExpressionR, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for a vector product multiplication
template<typename ExpressionL, typename ExpressionR, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionOuterProduct<ExpressionL, ExpressionR>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationOuterProduct <
    typename ExpressionToEvaluation<ExpressionL, testDim, solutionDim, spaceDim>::evaluation_Type
    , typename ExpressionToEvaluation<ExpressionR, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for a element-wise multiplication
template<typename ExpressionL, typename ExpressionR, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionEmult<ExpressionL, ExpressionR>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationEmult <
    typename ExpressionToEvaluation<ExpressionL, testDim, solutionDim, spaceDim>::evaluation_Type
    , typename ExpressionToEvaluation<ExpressionR, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for a division
template<typename ExpressionL, typename ExpressionR, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionDivision<ExpressionL, ExpressionR>, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationDivision <
    typename ExpressionToEvaluation<ExpressionL, testDim, solutionDim, spaceDim>::evaluation_Type
    , typename ExpressionToEvaluation<ExpressionR, testDim, solutionDim, spaceDim>::evaluation_Type
    > evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for h_K
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionHK, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationHK<spaceDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for meas_K
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionMeas, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationMeas<spaceDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionMeasBDCurrentFE, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationMeasBDCurrentFE<spaceDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for the position
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionPosition, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationPosition<spaceDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// Specialized for the normal
template<UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation<ExpressionNormal, testDim, solutionDim, spaceDim>
{
public:
    typedef EvaluationNormal<spaceDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

template<typename MeshType, typename MapType, UInt FESpaceDim, UInt testDim, UInt solutionDim, UInt spaceDim>
class ExpressionToEvaluation <
    ExpressionIfCrossed<MeshType, MapType, FESpaceDim>, testDim, solutionDim, spaceDim >
{
public:
    typedef EvaluationIfCrossed<MeshType, MapType, FESpaceDim> evaluation_Type;
private:
    ExpressionToEvaluation();
    ~ExpressionToEvaluation();
};

// \endcond

} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
