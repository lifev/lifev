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
 *   
 @file
 @brief This file contains the definition of the EvaluationScalarToVector class.
 
 @date 02/2014
 @author Davide Forti <davide.forti@epfl.ch>
 */

#ifndef EVALUTATION_SCALARTOVECTOR_HPP
#define EVALUTATION_SCALARTOVECTOR_HPP


#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/eta/expression/ExpressionScalarToVector.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>

namespace LifeV
{
    
namespace ExpressionAssembly
{

//! Evaluation for the transpose of another Evaluation
/*!
 @author Davide Forti <davide.forti@epfl.ch>
 
 */
    
template <typename EvaluationType, UInt FieldDim>
class EvaluationScalarToVector
{
public:
    
    
    //! @name Public Types
    //@{
    
    //! Type of the value returned by the 'operand' to be transposed
    //  typedef typename EvaluationType::return_Type return_Type;
    
    //! Type of the value returned by this class
    // NOTE: specialized for 3d problems
    typedef VectorSmall<FieldDim> return_Type;
    
    //@}
    
    
    //! @name Static constants
    //@{
    
    //! Flag for the global current FE
    const static flag_Type S_globalUpdateFlag;
    
    //! Flag for the test current FE
    const static flag_Type S_testUpdateFlag;
    
    //! Flag for the solution current FE
    const static flag_Type S_solutionUpdateFlag;
    
    //@}
    
    
    //! @name Constructors, destructor
    //@{
    
    //! Copy constructor
    EvaluationScalarToVector (const EvaluationScalarToVector& eval)
    : M_evaluation1 (eval.M_evaluation1),
    M_evaluation2 (eval.M_evaluation2),
    M_evaluation3 (eval.M_evaluation3),
    M_outputVector ( )
    {}
    
    //! Constructor from the corresponding expression
    template< typename Expression>
    explicit EvaluationScalarToVector (const ExpressionScalarToVector<Expression, FieldDim>& expression)
    : M_evaluation1 (expression.exprEx1() ),
    M_evaluation2 (expression.exprEx2() ),
    M_evaluation3 (expression.exprEx3() ),
    M_outputVector ( )
    {}
    
    //! Destructor
    ~EvaluationScalarToVector()
    {}
    
    //@}
    
    
    //! @name Methods
    //@{
    
    //! Internal update method
    void update (const UInt& iElement)
    {
        M_evaluation1.update (iElement);
        M_evaluation2.update (iElement);
        M_evaluation3.update (iElement);
    }
    
    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << " vector from scalar expression ";
        EvaluationType::display (out);
    }
    
    //@}
    
    
    //! @name Set Methods
    //@{
    
    //! Setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* globalCFE)
    {
        M_evaluation1.setGlobalCFE (globalCFE);
        M_evaluation2.setGlobalCFE (globalCFE);
        M_evaluation3.setGlobalCFE (globalCFE);
    }
    
    //! Setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* testCFE)
    {
        M_evaluation1.setTestCFE (testCFE);
        M_evaluation2.setTestCFE (testCFE);
        M_evaluation3.setTestCFE (testCFE);
    }
    
    //! Setter for the solution FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* solutionCFE)
    {
        M_evaluation1.setSolutionCFE (solutionCFE);
        M_evaluation2.setSolutionCFE (solutionCFE);
        M_evaluation3.setSolutionCFE (solutionCFE);
    }
    
    //! Setter for the quadrature rule
    void setQuadrature (const QuadratureRule& qr)
    {
        M_evaluation1.setQuadrature (qr);
        M_evaluation2.setQuadrature (qr);
        M_evaluation3.setQuadrature (qr);
    }
    
    //@}
    
    
    //! @name Get Methods
    //@{
    
    //! Getter for a value
    return_Type value_q (const UInt& q) const
    {
        M_outputVector[0]  = M_evaluation1.value_q( q ) ;
        M_outputVector[1]  = M_evaluation2.value_q( q ) ;
        M_outputVector[2]  = M_evaluation3.value_q( q ) ;
        return M_outputVector;
    }
    
    //! Getter for the value for a vector
    return_Type value_qi (const UInt& q, const UInt& i) const
    {
        
        M_outputVector[0] = M_evaluation1.value_q( q ,i) ;
        M_outputVector[1] = M_evaluation2.value_q( q ,i) ;
        M_outputVector[2] = M_evaluation3.value_q( q ,i) ;
        return M_outputVector;
    }
    
    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& q, const UInt& i, const UInt& j) const
    {
        M_outputVector[0] = M_evaluation1.value_q( q,i,j ) ;
        M_outputVector[1] = M_evaluation2.value_q( q,i,j ) ;
        M_outputVector[2] = M_evaluation3.value_q( q,i,j ) ;
        return M_outputVector;
    }
    
    //@}
    
private:
    
    //! @name Private Methods
    //@{
    
    //! No default
    EvaluationScalarToVector();
    
    //@}
    
    // Internal storage
    EvaluationType M_evaluation1;
    EvaluationType M_evaluation2;
    EvaluationType M_evaluation3;
    VectorSmall<FieldDim> M_outputVector;
};



template< typename EvaluationType, UInt FieldDim >
const flag_Type EvaluationScalarToVector<EvaluationType, FieldDim>::S_globalUpdateFlag
= EvaluationType::S_globalUpdateFlag;

template< typename EvaluationType, UInt FieldDim >
const flag_Type EvaluationScalarToVector<EvaluationType, FieldDim>::S_testUpdateFlag
= EvaluationType::S_testUpdateFlag;

template< typename EvaluationType, UInt FieldDim >
const flag_Type EvaluationScalarToVector<EvaluationType, FieldDim>::S_solutionUpdateFlag
= EvaluationType::S_solutionUpdateFlag;

} // Namespace ExpressionAssembly
    
} // Namespace LifeV
#endif
