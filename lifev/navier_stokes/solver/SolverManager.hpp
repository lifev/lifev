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

 */

#ifndef SOLVERMANAGER_HPP_
#define SOLVERMANAGER_HPP_

#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <lifev/core/util/Factory.hpp>

#include <lifev/core/linear_algebra/InvertibleOperator.hpp>
#include <lifev/core/linear_algebra/BlockOperator.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

#include <Teuchos_ParameterList.hpp>

#include <lifev/core/util/LifeChrono.hpp>

#include <lifev/navier_stokes/solver/aSIMPLE.hpp>

namespace LifeV
{

class SolverManager
{
public:

    //! @name Public Typedefs
    //@{

    typedef Operators::LinearOperator::comm_Type    comm_Type;
    typedef Operators::LinearOperator::commPtr_Type commPtr_Type;

    typedef Operators::LinearOperator::operator_Type operator_Type;
    typedef Operators::LinearOperator::operatorPtr_Type operatorPtr_Type;

    typedef Operators::LinearOperator::vector_Type  vector_Type;

    typedef Teuchos::ParameterList parameterList_Type;
    typedef boost::shared_ptr<parameterList_Type> parameterListPtr_Type;
    
    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type > matrixPtr_Type;
    
    typedef Epetra_CrsMatrix rowMatrix_Type;
    typedef boost::shared_ptr<rowMatrix_Type> rowMatrixPtr_Type;

    //@}

    SolverManager();

    virtual ~SolverManager(){};

    void setComm(const commPtr_Type & comm);

    void setSchurOptions(const parameterListPtr_Type & schurComplementList);

    void setMomentumOptions(const parameterListPtr_Type & momentumList);

    void setLinSolverParameter(const parameterListPtr_Type & linearSolverList);

    const operatorPtr_Type updateInvertibleOperator( const matrixPtr_Type& F, const matrixPtr_Type& B, const matrixPtr_Type& Btranspose );

protected:

    virtual void getMatrices( const matrixPtr_Type& F, const matrixPtr_Type& B, const matrixPtr_Type& Btranspose ) = 0;
    
    virtual void updateApproximatedMomentumOperator( ) = 0;

    virtual void updateApproximatedSchurComplementOperator( ) = 0;

    rowMatrixPtr_Type M_F;
    rowMatrixPtr_Type M_Btranspose;
    rowMatrixPtr_Type M_B;
    
    operatorPtr_Type M_approximatedSchurComplementOperator;
    operatorPtr_Type M_approximatedMomentumOperator;

    boost::shared_ptr<Operators::BlockOperator> M_oper;
    boost::shared_ptr<Operators::BlockOperator> M_prec;
    boost::shared_ptr<Operators::InvertibleOperator> M_invOper;

    parameterListPtr_Type M_linearSolverList;
    parameterListPtr_Type M_momentumList;
    parameterListPtr_Type M_schurComplementList;

    Displayer M_displayer;

};

typedef FactorySingleton<Factory<SolverManager, std::string> > SolverFactory;

} /*end namespace */


#endif /* SOLVERMANAGER_HPP_ */
