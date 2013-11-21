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
    @file SolverPolicyLinearSolver class
    @brief This class contains all the informations necessary
           to solve a linear system

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 07-12-2012
 */

#include <lifev/navier_stokes/solver/NavierStokesSolver/SolverPolicyLinearSolver.hpp>

#include <string>
#include <lifev/core/util/LifeChrono.hpp>

namespace LifeV
{

void
SolverPolicyLinearSolver::initSolver ( Teuchos::ParameterList& list )
{
    // Parameter list for the solver
    std::string resourcesPath         = list.get ( "Resources path", "./Resources" );
    resourcesPath.append ("/");
    const std::string solverParamFile = list.get ( "Parameters file", "SolverParamList.xml" );
    Teuchos::RCP< Teuchos::ParameterList > solverParameters = Teuchos::rcp ( new Teuchos::ParameterList );
    solverParameters = Teuchos::getParametersFromXmlFile ( resourcesPath + solverParamFile );

    // Solver initialization
    displayer().leaderPrint ( "\n[Solver initialization]\n" );
    M_solver.reset ( new solver_Type );
    M_solver->setCommunicator ( comm() );

    displayer().leaderPrint ( "Setting up the solver... " );
    M_solver->setParameters ( *solverParameters );
    M_solver->showMe();
}

void
SolverPolicyLinearSolver::setPreconditioner ( preconditionerPtr_Type preconditionerPtr )
{
    M_solver->setPreconditioner ( preconditionerPtr );
}

SolverPolicyLinearSolver::preconditionerPtr_Type
SolverPolicyLinearSolver::preconditioner()
{
    return M_solver->preconditioner();
}


int
SolverPolicyLinearSolver::solve ( matrixPtr_Type systemMatrix,
                                  vectorPtr_Type rhs,
                                  vectorPtr_Type solution )
{
    M_solver->setOperator ( systemMatrix );
    M_solver->setRightHandSide ( rhs );
    return M_solver->solve ( solution );
}

} // namespace LifeV
