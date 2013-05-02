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
    @brief AztecooOperator

    @author Umberto Villa <umberto.villa@gmail.com>

    @date 03-09-2010
 */

#ifndef _AZTECOOOPERATOR_HPP_
#define _AZTECOOOPERATOR_HPP_

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"

#include <AztecOO.h>
#include <Teuchos_ParameterList.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wextra"

#include <lifev/core/operator/SolverOperator.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{
namespace Operators
{
//! @class AztecooOperator
/*! @brief Abstract class which defines the interface of an Invertible Linear Operator.
 *
 */
class AztecooOperator : public SolverOperator
{
public:
    typedef AztecOO SolverType;
    typedef boost::shared_ptr<SolverType> SolverType_ptr;

    AztecooOperator();

protected:

    virtual int doApplyInverse ( const vector_Type& X, vector_Type& Y ) const;
    virtual void doSetOperator();
    virtual void doSetPreconditioner();
    virtual void doSetParameterList();
    virtual void doResetSolver();

    SolverType_ptr M_linSolver;
};

inline SolverOperator* createAztecooOperator()
{
    return new AztecooOperator();
}
namespace
{
static bool registerAztecoo = SolverOperatorFactory::instance().registerProduct ( "AztecOO", &createAztecooOperator );
}


} // Namespace Operators

} // Namespace LifeV

#endif // _AZTECOOOPERATOR_HPP_
