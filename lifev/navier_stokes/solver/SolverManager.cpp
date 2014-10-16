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

#include <lifev/core/util/LifeChrono.hpp>

#include <lifev/navier_stokes/solver/SolverManager.hpp>

#include <lifev/core/linear_algebra/LumpedOperator.hpp>

#include <lifev/core/linear_algebra/IfpackPreconditioner.hpp>
#include <lifev/core/linear_algebra/MLPreconditioner.hpp>
#include <lifev/core/linear_algebra/TwoLevelPreconditioner.hpp>
#include <lifev/core/linear_algebra/AztecooOperator.hpp>
#include <lifev/core/linear_algebra/BelosOperator.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>

#include <ml_epetra_utils.h>

namespace LifeV
{

SolverManager::SolverManager():
    M_oper(new Operators::BlockOperator),
    M_prec(new Operators::BlockOperator),
    M_invOper(),
    M_displayer()
{
}

//===========================================================================//
//                          Setters methods                                  //
//===========================================================================//

void SolverManager::setComm(const commPtr_Type & comm)
{
    M_displayer.setCommunicator(comm);
}


void SolverManager::setMomentumOptions(const parameterListPtr_Type & _oList)
{
    //M_momentumOptions = _oList;
}


void SolverManager::setSchurOptions(const parameterListPtr_Type & _oList)
{
    //M_schurOptions = _oList;
}

void SolverManager::setLinSolverParameter(const parameterListPtr_Type & _pList)
{
    //M_pListLinSolver = _pList;
}


//const SolverManager::operatorPtr_Type SolverManager::updateInvertibleOperator( )
//{
	/*
    LifeChrono chrono;

    //(1) Do Manager specific operations
    M_displayer.leaderPrint( "OseenOperatorManager:\n\tparsing matrix container...");
    chrono.start();
    parseMatrixContainer( _mc );
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    M_displayer.leaderPrint( "\tupdate the approximate momentum operator...");
    chrono.start();
    updateApproximatedMomentumOperator( );
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    M_displayer.leaderPrint( "\tupdate the approximate Schur complement operator...");
    chrono.start();
    updateApproximatedSchurComplementOperator( );
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    ASSERT_POS(M_momentumMatrix.get() != 0,
               "[OseenOperatorManager::updateInvertibleOperator]: parseMatrixContainer did not initialize M_momentumMatrix");
    ASSERT_POS(M_pressureGradientMatrix.get() != 0,
               "[OseenOperatorManager::updateInvertibleOperator]: parseMatrixContainer did not initialize M_pressureGradientMatrix");
    ASSERT_POS(M_velocityDivergenceMatrix.get() != 0,
               "[OseenOperatorManager::updateInvertibleOperator]: parseMatrixContainer did not initialize M_velocityDivergenceMatrix");
    ASSERT_POS(M_approximatedMomentumOperator.get() != 0,
               "[OseenOperatorManager::updateInvertibleOperator]: updateApproximatedMomentumOperator did not initialize M_approximatedMomentumOperator");
    ASSERT_POS(M_approximatedSchurComplementOperator.get() != 0,
               "[OseenOperatorManager::updateInvertibleOperator]:updateApproximatedSchurComplementOperator did not initialize M_approximatedSchurComplementOperator");

    //(2) Set up the OseenOperator
    M_displayer.leaderPrint( "\tset up the block operator...");
    chrono.start();

    Operators::BlockOperator::operatorPtrContainer_Type operData(2,2);
    operData(0,0) = M_momentumMatrix;
    operData(0,1) = M_pressureGradientMatrix;
    operData(1,0) = M_velocityDivergenceMatrix;

    M_oper->setUp(operData, M_displayer.comm());
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    //(3) Set up the OseenPreconditioner
    M_displayer.leaderPrint( "\tset up the block preconditioner...");
    chrono.start();
    Operators::BlockOperator::operatorPtrContainer_Type precData(2,2);
    precData(0,0) = M_approximatedMomentumOperator;
    precData(0,1) = M_pressureGradientMatrix;
    precData(1,1) = M_approximatedSchurComplementOperator;
    M_prec->setUp(precData, M_displayer.comm());
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    //(4) Set up the OseenSolver
    M_displayer.leaderPrint( "\tset up the Trilinos solver...");
    chrono.start();
    std::string solverType(M_pListLinSolver->get<std::string>("Linear Solver Type"));
    M_invOper.reset(Operators::InvertibleOperatorFactory::instance().createObject(solverType));

    M_invOper->setOperator(M_oper);
    M_invOper->setPreconditioner(M_prec);
    M_invOper->setParameterList(M_pListLinSolver->sublist(solverType));

    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );
	*/

//    return M_invOper;
//}

} /*end namespace */
