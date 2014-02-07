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
 @brief This expression will collect three expressions that return a scalar into a vector. Used for computing the laplacian of a vector field.
 
 @author Davide Forti <davide.forti@epfl.ch>
 
 @date 02-2014
 */

#ifndef EXPRESSION_SCALARTOVECTOR_HPP
#define EXPRESSION_SCALARTOVECTOR_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>

namespace LifeV
{
    
namespace ExpressionAssembly
{

template <typename ExpressionType, UInt FieldDim>
class ExpressionScalarToVector : public ExpressionBase<ExpressionScalarToVector<ExpressionType, FieldDim> >
{
public:
    
    //! @name Public Types
    //@{
    
    // No real need, just for ease of coding
    typedef ExpressionBase< ExpressionScalarToVector <ExpressionType, FieldDim> > base_Type;
    
    //@}
    
    
    //! @name Constructors & Destructor
    //@{
    
    //! Full constructor
    ExpressionScalarToVector(const ExpressionType& expr1, const ExpressionType& expr2, const ExpressionType& expr3)
    : base_Type(), M_expr1 (expr1), M_expr2 (expr2), M_expr3 (expr3) {}
    
    //! Copy constructor
    ExpressionScalarToVector (const ExpressionScalarToVector<ExpressionType, FieldDim>& expression)
    : base_Type(), M_expr1 (expression.M_expr1), M_expr2 (expression.M_expr2), M_expr3 (expression.M_expr3)  {}
    
    //! Destructor
    ~ExpressionScalarToVector() {}
    
    //@}
    
    
    //! @name Methods
    //@{
    
    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << " vector small from expression ";
        ExpressionType::display (out);
    }
    
    //@}
    
    
    //! @name Get Methods
    //@{
    
    //! Getter for the expression that we transpose
    const ExpressionType& exprEx1() const
    {
        return M_expr1;
    }
    
    const ExpressionType& exprEx2() const
    {
        return M_expr2;
    }
    
    const ExpressionType& exprEx3() const
    {
        return M_expr3;
    }
    
    //@}
    
private:
    
    //! @name Private Methods
    //@{
    
    //! No default constructor
    ExpressionScalarToVector();
    
    //@}
    
    // Expression that we want to collect into the VectorSmall
    ExpressionType M_expr1;
    ExpressionType M_expr2;
    ExpressionType M_expr3;
    
    
};

template< typename ExpressionType, UInt FieldDim >
ExpressionScalarToVector<ExpressionType, FieldDim>
scalarToVector(const ExpressionBase<ExpressionType>& expr1,
                         const ExpressionBase<ExpressionType>& expr2,
                         const ExpressionBase<ExpressionType>& expr3)
{
    return ExpressionScalarToVector<ExpressionType, FieldDim> (expr1.cast(), expr2.cast(), expr3.cast());
}
    
    
} // Namespace ExpressionAssembly
    
} // Namespace LifeV
#endif