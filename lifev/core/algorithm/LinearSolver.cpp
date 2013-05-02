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
    @brief LinearSolver

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 03-08-2011
 */

#include <lifev/core/LifeV.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
LinearSolver::LinearSolver() :
    M_operator             (),
    M_matrix               (),
    M_baseMatrixForPreconditioner(),
    M_preconditioner       (),
    M_preconditionerOperator(),
    M_solverType           ( UndefinedSolver ),
    M_solverOperator       (),
    M_parameterList        (),
    M_displayer            ( new Displayer() ),
    M_maxItersForReuse     ( 0 ),
    M_reusePreconditioner  ( false ),
    M_quitOnFailure        ( false ),
    M_silent               ( false ),
    M_lossOfPrecision      ( SolverOperator_Type::undefined ),
    M_maxNumItersReached   ( SolverOperator_Type::undefined ),
    M_converged            ( SolverOperator_Type::undefined ),
    M_tolerance            ( -1. )
{

}

LinearSolver::LinearSolver ( const boost::shared_ptr<Epetra_Comm> commPtr ) :
    M_operator             (),
    M_matrix               (),
    M_baseMatrixForPreconditioner(),
    M_preconditioner       (),
    M_preconditionerOperator(),
    M_solverType           ( UndefinedSolver ),
    M_solverOperator       (),
    M_parameterList        (),
    M_displayer            ( new Displayer ( commPtr ) ),
    M_maxItersForReuse     ( 0 ),
    M_reusePreconditioner  ( false ),
    M_quitOnFailure        ( false ),
    M_silent               ( false ),
    M_lossOfPrecision      ( SolverOperator_Type::undefined ),
    M_maxNumItersReached   ( SolverOperator_Type::undefined ),
    M_converged            ( SolverOperator_Type::undefined ),
    M_tolerance            ( -1. )
{

}

LinearSolver::~LinearSolver()
{

}

// ===================================================
// Methods
// ===================================================
Int
LinearSolver::solve ( vectorPtr_Type solutionPtr )
{
    // Build preconditioners if needed
    bool retry ( true );
    if ( !isPreconditionerSet() || !M_reusePreconditioner  )
    {
        buildPreconditioner();

        // There will be no retry if the preconditioner is recomputed
        retry = false;
    }
    else
    {
        if ( !M_silent )
        {
            M_displayer->leaderPrint ( "SLV-  Reusing precond ...\n" );
        }
    }

    if ( M_rhs.get() == 0 || M_operator == 0 )
    {
        M_displayer->leaderPrint ( "SLV-  ERROR: LinearSolver failed to set up correctly!\n" );
        return -1;
    }

    // Setup the Solver Operator
    setupSolverOperator();

    // Reset status informations
    bool failure = false;
    this->resetStatus();
    M_solverOperator->resetStatus();

    // Solve the linear system
    LifeChrono chrono;
    chrono.start();

    M_solverOperator->ApplyInverse ( M_rhs->epetraVector(), solutionPtr->epetraVector() );
    M_converged         = M_solverOperator->hasConverged();
    M_lossOfPrecision   = M_solverOperator->isLossOfAccuracyDetected();
    chrono.stop();
    if ( !M_silent )
    {
        M_displayer->leaderPrintMax ( "SLV-  Solution time: " , chrono.diff(), " s." );
    }

    // Getting informations post-solve
    Int numIters = M_solverOperator->numIterations();

    // Second run recomputing the preconditioner
    // This is done only if the preconditioner has not been
    // already recomputed and if it is a LifeV preconditioner.
    if ( M_converged != SolverOperator_Type::yes
            && retry
            && M_preconditioner )
    {
        M_displayer->leaderPrint ( "SLV-  Iterative solver failed, numiter = " , numIters, "\n" );
        M_displayer->leaderPrint ( "SLV-  retrying:\n" );

        buildPreconditioner();

        // Solving again, but only once (retry = false)
        chrono.start();
        M_solverOperator->ApplyInverse ( M_rhs->epetraVector(), solutionPtr->epetraVector() );
        M_converged         = M_solverOperator->hasConverged();
        M_lossOfPrecision   = M_solverOperator->isLossOfAccuracyDetected();
        chrono.stop();
        if ( !M_silent )
        {
            M_displayer->leaderPrintMax ( "SLV-  Solution time: " , chrono.diff(), " s." );
        }
    }

    if ( M_lossOfPrecision == SolverOperator_Type::yes )
    {
        M_displayer->leaderPrint ( "SLV-  WARNING: Loss of accuracy detected!\n" );
        failure = true;
    }

    if ( M_converged == SolverOperator_Type::yes )
    {
        if ( !M_silent )
        {
            M_displayer->leaderPrint ( "SLV-  Convergence in " , numIters, " iterations\n" );
        }
        M_maxNumItersReached = SolverOperator_Type::no;
    }
    else
    {
        M_displayer->leaderPrint ( "SLV-  WARNING: Solver failed to converged to the desired precision!\n" );
        M_maxNumItersReached = SolverOperator_Type::yes;
        failure = true;
    }

    // If quitOnFailure is enabled and if some problems occur
    // the simulation is stopped
    if ( M_quitOnFailure && failure )
    {
        exit ( -1 );
    }

    // Reset the solver to free the internal pointers
    M_solverOperator->resetSolver();

    // If the number of iterations reaches the threshold of maxIterForReuse
    // we reset the preconditioners to force to solver to recompute it next
    // time
    if ( numIters > M_maxItersForReuse )
    {
        resetPreconditioner();
    }

    return numIters;
}

Real
LinearSolver::computeResidual ( vectorPtr_Type solutionPtr )
{
    if ( !M_operator || !M_rhs )
    {
        M_displayer->leaderPrint ( "SLV-  WARNING: LinearSolver can not compute the residual if the operator and the RHS are not set!\n" );
        return -1;
    }

    vector_Type Ax ( solutionPtr->map() );
    vector_Type residual ( *M_rhs );

    M_operator->Apply ( solutionPtr->epetraVector(), Ax.epetraVector() );

    residual.epetraVector().Update ( 1, Ax.epetraVector(), -1 );

    Real residualNorm;

    residual.norm2 ( &residualNorm );

    return residualNorm;
}

std::string
LinearSolver::printStatus()
{
    std::ostringstream stat;
    std::string str;

    if ( M_lossOfPrecision == SolverOperator_Type::yes )
    {
        stat << "Accuracy loss ";
    }
    if ( M_maxNumItersReached == SolverOperator_Type::yes )
    {
        stat << "Maximum number of iterations reached ";
    }
    if ( M_converged == SolverOperator_Type::yes )
    {
        stat << "The solver has converged ";
    }
    else if ( M_converged == SolverOperator_Type::no )
    {
        stat << "The solver has not ";
    }
    str = stat.str();
    return str;
}

void
LinearSolver::setPreconditionerFromGetPot ( const GetPot& dataFile, const std::string& section )
{
    std::string precName = dataFile ( ( section + "/prectype" ).data(), "Ifpack" );

    M_preconditioner.reset ( PRECFactory::instance().createObject ( precName ) );
    ASSERT ( M_preconditioner.get() != 0, " Preconditioner not set" );

    M_preconditioner->setDataFromGetPot ( dataFile, section );
}

void
LinearSolver::buildPreconditioner()
{
    LifeChrono chrono;
    Real condest ( -1 );

    if ( M_preconditioner )
    {
        if ( M_matrix.get() == 0 )
        {
            M_displayer->leaderPrint ( "SLV-  ERROR: LinearSolver requires a matrix to build the preconditioner!\n" );
            exit ( 1 );
        }
        else
        {
            chrono.start();
            if ( !M_silent )
            {
                M_displayer->leaderPrint ( "SLV-  Computing the preconditioner...\n" );
            }
            if ( M_baseMatrixForPreconditioner.get() == 0 )
            {
                if ( !M_silent )
                {
                    M_displayer->leaderPrint ( "SLV-  Build the preconditioner using the problem matrix\n" );
                }
                M_preconditioner->buildPreconditioner ( M_matrix );
            }
            else
            {
                if ( !M_silent )
                {
                    M_displayer->leaderPrint ( "SLV-  Build the preconditioner using the base matrix provided\n" );
                }
                M_preconditioner->buildPreconditioner ( M_baseMatrixForPreconditioner );
            }
            condest = M_preconditioner->condest();
            chrono.stop();
            if ( !M_silent )
            {
                M_displayer->leaderPrintMax ( "SLV-  Preconditioner computed in " , chrono.diff(), " s." );
            }
            if ( !M_silent )
            {
                M_displayer->leaderPrint ( "SLV-  Estimated condition number               " , condest, "\n" );
            }
        }
    }
}

void
LinearSolver::resetPreconditioner()
{
    if ( M_preconditioner )
    {
        M_preconditioner->resetPreconditioner();
    }
}

bool
LinearSolver::isPreconditionerSet() const
{
    if ( M_preconditionerOperator )
    {
        return true;
    }

    return M_preconditioner.get() != 0 && M_preconditioner->preconditionerCreated();
}

void
LinearSolver::resetStatus()
{
    M_lossOfPrecision    = SolverOperator_Type::undefined;
    M_maxNumItersReached = SolverOperator_Type::undefined;
    M_converged          = SolverOperator_Type::undefined;
}

void
LinearSolver::showMe ( std::ostream& output ) const
{
    if ( M_displayer->isLeader() )
    {
        output << "Solver parameters list:" << std::endl;
        output << "-----------------------------" << std::endl;
        output << M_parameterList << std::endl;
        output << "-----------------------------" << std::endl;
    }
}

void
LinearSolver::setupSolverOperator()
{
    // Creation of a solver if there exists any
    if ( !M_solverOperator )
    {
        switch ( M_solverType )
        {
            case Belos:
                M_solverOperator.reset ( Operators::SolverOperatorFactory::instance().createObject ( "Belos" ) );
                break;
            case AztecOO:
                M_solverOperator.reset ( Operators::SolverOperatorFactory::instance().createObject ( "AztecOO" ) );
                break;
            default:
                M_displayer->leaderPrint ( "SLV-  ERROR: The type of solver is not recognized!\n" );
                exit ( 1 );
                break;
        }
    }

    // Set the preconditioner operator in the SolverOperator object
    if ( M_preconditioner )
    {
        M_solverOperator->setPreconditioner ( M_preconditioner->preconditionerPtr() );
    }
    else if ( M_preconditionerOperator )
    {
        M_solverOperator->setPreconditioner ( M_preconditionerOperator );
    }

    // Set the operator in the SolverOperator object
    M_solverOperator->setOperator ( M_operator );

    // Set the tolerance if it has been set
    if ( M_tolerance > 0 )
    {
        M_solverOperator->setTolerance ( M_tolerance );
    }

    // Set the parameter inside the solver
    M_solverOperator->setParameters ( M_parameterList.sublist ( "Solver: Operator List" ) );
}

// ===================================================
// Set Methods
// ===================================================
void
LinearSolver::setSolverType ( const SolverType& solverType )
{
    M_solverType = solverType;
}

void
LinearSolver::setCommunicator ( const boost::shared_ptr<Epetra_Comm> commPtr )
{
    M_displayer->setCommunicator ( commPtr );
}

void LinearSolver::setOperator ( matrixPtr_Type matrixPtr )
{
    M_operator = matrixPtr->matrixPtr();
    M_matrix = matrixPtr;
}

void
LinearSolver::setOperator ( operatorPtr_Type operPtr )
{
    M_matrix.reset();
    M_operator = operPtr;
}

void
LinearSolver::setRightHandSide ( const vectorPtr_Type rhsPtr )
{
    M_rhs = rhsPtr;
}

void
LinearSolver::setPreconditioner ( preconditionerPtr_Type preconditionerPtr )
{
    // If a preconditioner operator exists it must be deleted
    M_preconditionerOperator.reset();

    M_preconditioner = preconditionerPtr;
}

void
LinearSolver::setPreconditioner ( operatorPtr_Type preconditionerPtr )
{
    // If a LifeV::Preconditioner exists it must be deleted
    M_preconditioner.reset();

    M_preconditionerOperator = preconditionerPtr;
}

void
LinearSolver::setBaseMatrixForPreconditioner ( matrixPtr_Type baseMatrixPtr )
{
    M_baseMatrixForPreconditioner = baseMatrixPtr;
}

void
LinearSolver::setParameters ( const Teuchos::ParameterList& list )
{
    M_parameterList.setParameters ( list );
    M_solverType = UndefinedSolver;
    std::string solverName = M_parameterList.get<std::string> ( "Solver Type" );
    if ( solverName == "Belos" )
    {
        M_solverType = Belos;
    }
    else if ( solverName == "AztecOO" )
    {
        M_solverType = AztecOO;
    }

    M_reusePreconditioner  = M_parameterList.get ( "Reuse Preconditioner"     , false );
    Int maxIter            = M_parameterList.get ( "Maximum Iterations"       , 200 );
    M_maxItersForReuse     = M_parameterList.get ( "Max Iterations For Reuse" , static_cast<Int> ( maxIter * 8. / 10. ) );
    M_quitOnFailure        = M_parameterList.get ( "Quit On Failure"          , false );
    M_silent               = M_parameterList.get ( "Silent"                   , false );
}

void
LinearSolver::resetParameters()
{
    M_parameterList = Teuchos::ParameterList();
}

void
LinearSolver::setReusePreconditioner ( const bool reusePreconditioner )
{
    M_reusePreconditioner = reusePreconditioner;
}

void
LinearSolver::setQuitOnFailure ( const bool enable )
{
    M_quitOnFailure = enable;
}

void
LinearSolver::setTolerance ( const Real& tolerance )
{
    M_tolerance = tolerance;
}

// ===================================================
// Get Methods
// ===================================================

Int
LinearSolver::numIterations() const
{
    return M_solverOperator->numIterations();
}


Real
LinearSolver::recursiveResidual()
{
    M_displayer->leaderPrint ( "SLV-  WARNING: LinearSoler::recursiveResidual is not yet implemented\n" );

    /*
    if ( !M_problem->isProblemSet() )
    {
        M_displayer->leaderPrint( "SLV-  WARNING: LinearSoler can not compute the residual if the linear system is not set!\n" );
        return -1;
    }
    multiVector_Type res( *( M_problem->getRHS().get() ) );
    M_problem->computeCurrResVec( &res,M_problem->getLHS().get(), M_problem->getRHS().get() );
    Real residual;
    res.Norm2( &residual );
    return residual;
    */
    return 0.;
}

LinearSolver::preconditionerPtr_Type&
LinearSolver::preconditioner()
{
    return M_preconditioner;
}

Teuchos::ParameterList&
LinearSolver::parametersList()
{
    return M_parameterList;
}

LinearSolver::SolverOperatorPtr_Type
LinearSolver::solver()
{
    return M_solverOperator;
}

boost::shared_ptr<Displayer>
LinearSolver::displayer()
{
    return M_displayer;
}

Int
LinearSolver::maxItersForReuse() const
{
    return M_maxItersForReuse;
}

bool
LinearSolver::reusePreconditioner() const
{
    return M_reusePreconditioner;
}

bool
LinearSolver::quitOnFailure() const
{
    return M_quitOnFailure;
}

bool
LinearSolver::silent() const
{
    return M_silent;
}

LinearSolver::SolverOperator_Type::SolverOperatorStatusType
LinearSolver::hasReachedMaxNumIters() const
{
    return M_maxNumItersReached;
}

LinearSolver::SolverOperator_Type::SolverOperatorStatusType
LinearSolver::isLossOfAccuracyDetected() const
{
    return M_lossOfPrecision;
}

LinearSolver::SolverOperator_Type::SolverOperatorStatusType
LinearSolver::hasConverged() const
{
    return M_converged;
}

// ===================================================
// Private Methods
// ===================================================

} // namespace LifeV
