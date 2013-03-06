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
    @brief BelosOperator

    @author Umberto Villa <umberto.villa@gmail.com>

    @date 28-09-2010
 */

#ifndef _BELOSOPERATOR_HPP_
#define _BELOSOPERATOR_HPP_

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"

#include <BelosEpetraAdapter.hpp>
#include <BelosSolverManager.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RefCountPtr.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wextra"

#include <lifev/core/operator/SolverOperator.hpp>

namespace LifeV
{
namespace Operators
{
//! @class BelosOperator
/*! @brief Class which defines the interface of an Invertible Linear Operator through belos.
 *
 */

class BelosOperator : public SolverOperator
{
public:

    //! @name Public Typedefs and Enumerators
    //@{

    enum PreconditionerSide { None, Left, Right };

    enum SolverManagerType { NotAValidSolverManager, BlockCG, PseudoBlockCG, RCG,
                             BlockGmres, PseudoBlockGmres, GmresPoly,
                             GCRODR, PCPG, TFQMR, MINRES
                           };

    //@}

    //! null constructor and destructor
    //@{
    BelosOperator();
    ~BelosOperator();
    //@}

protected:

    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator    OP;
    typedef Belos::LinearProblem<double, MV, OP> LinearProblem;
    typedef Belos::SolverManager<double, MV, OP> SolverType;
    typedef Teuchos::RCP<LinearProblem> LinearProblem_ptr;
    typedef Teuchos::RCP<SolverType>    SolverType_ptr;


    virtual int doApplyInverse ( const vector_Type& X, vector_Type& Y ) const;
    virtual void doSetOperator();
    virtual void doSetPreconditioner();
    virtual void doSetParameterList();
    virtual void doResetSolver();
    void allocateSolver ( const SolverManagerType& solverManagerType );
    //! The linearProblem
    LinearProblem_ptr M_linProblem;
    //! The linearSolver
    SolverType_ptr M_solverManager;
    //! Cast to a Belos Preconditioner
    Teuchos::RCP<Belos::EpetraPrecOp> M_belosPrec;

    static SolverManagerType  getSolverManagerTypeFromString ( const std::string& str );
    static PreconditionerSide getPreconditionerSideFromString ( const std::string& str );

};

inline SolverOperator* createBelosOperator()
{
    return new BelosOperator();
}
namespace
{
static bool registerBelos = SolverOperatorFactory::instance().registerProduct ( "Belos", &createBelosOperator );
}


} /*end namespace Operators */
} /*end namespace LifeV */
#endif /* _BELOSOPERATOR_HPP_ */
