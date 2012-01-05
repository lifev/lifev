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

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include "Teuchos_RCPBoostSharedPtrConversions.hpp"

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/LifeV.hpp>
#include <life/lifealg/LinearSolver.hpp>
#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
LinearSolver::LinearSolver() :
        M_leftPreconditioner   (),
        M_rightPreconditioner  (),
        M_solverManagerType    ( BlockGmres ),
        M_solverManager        (),
        M_problem              ( new LinearProblem_Type() ),
        M_parameterList        (),
        M_displayer            ( new Displayer() ),
        M_maxItersForReuse     ( 0 ),
        M_reusePreconditioner  ( false ),
        M_quitOnFailure        ( false ),
        M_silent               ( false ),
        M_lossOfPrecision      ( false ),
        M_maxNumItersReached   ( false )
{
    M_problem->setLabel( "SolverBelos" );
}

LinearSolver::LinearSolver( const boost::shared_ptr<Epetra_Comm>& comm ) :
        M_leftPreconditioner   (),
        M_rightPreconditioner  (),
        M_solverManagerType    ( BlockGmres ),
        M_solverManager        (),
        M_problem              ( new LinearProblem_Type() ),
        M_parameterList        (),
        M_displayer            ( new Displayer( comm ) ),
        M_maxItersForReuse     ( 0 ),
        M_reusePreconditioner  ( false ),
        M_quitOnFailure        ( false ),
        M_silent               ( false ),
        M_lossOfPrecision      ( false ),
        M_maxNumItersReached   ( false )
{
    M_problem->setLabel( "SolverBelos" );
}

LinearSolver::~LinearSolver()
{

}

// ===================================================
// Methods
// ===================================================
Int
LinearSolver::solve( vector_Type& solution )
{
    return this->solve( solution.epetraVector() );
}

Int
LinearSolver::solve( multiVector_Type& solution )
{
    // Reset status informations
    M_lossOfPrecision = false;
    M_maxNumItersReached = false;
    bool failure = false;

    // Setting the unknown in the system
    Teuchos::RCP<multiVector_Type> solutionPtr( &( solution ), false );
    M_problem->setLHS( solutionPtr );

    // Build preconditioners if needed
    bool retry( true );
    if ( !isPreconditionerSet() || !M_reusePreconditioner  )
    {
        buildPreconditioner();

        // There will be no retry if the preconditioner is recomputed
        retry = false;
    }
    else
    {
        if( !M_silent ) M_displayer->leaderPrint( "SLV-  Reusing precond ...\n" );
    }

    bool set = M_problem->setProblem();
    if ( set == false ) {
        M_displayer->leaderPrint( "SLV-  ERROR: SolverBelos failed to set up correctly!\n" );
        return -1;
    }

    // Setup the Solver Manager
    setupSolverManager();

    // Solve the linear system
    LifeChrono chrono;
    chrono.start();
    Belos::ReturnType ret = M_solverManager->solve();
    chrono.stop();
    if( !M_silent ) M_displayer->leaderPrintMax( "SLV-  Solution time: " , chrono.diff(), " s." );

    // Getting informations post-solve
    Int numIters = M_solverManager->getNumIters();
    bool M_lossOfPrecision = M_solverManager->isLOADetected();

    // Second run recomputing the preconditioner
    // This is done only if the preconditioner has not been
    // already recomputed and if it is a LifeV preconditioner.
    if ( ret != Belos::Converged && retry &&
       ( M_leftPreconditioner || M_rightPreconditioner ) )
    {
        M_displayer->leaderPrint( "SLV-  Iterative solver failed, numiter = " , numIters, "\n" );
        M_displayer->leaderPrint( "SLV-  retrying:\n" );

        buildPreconditioner();

        // Solving again, but only once (retry = false)
        chrono.start();
        ret = M_solverManager->solve();
        chrono.stop();
        if( !M_silent ) M_displayer->leaderPrintMax( "SLV-  Solution time: " , chrono.diff(), " s." );
    }

    if ( M_lossOfPrecision )
    {
        M_displayer->leaderPrint( "SLV-  WARNING: Loss of accuracy detected!\n" );
        failure = true;
    }

    if ( ret == Belos::Converged )
    {
        if( !M_silent ) M_displayer->leaderPrint( "SLV-  Convergence in " , numIters, " iterations\n" );
    }
    else
    {
        M_displayer->leaderPrint( "SLV-  WARNING: SolverBelos failed to converged to the desired precision!\n" );
        M_maxNumItersReached = true;
        failure = true;
    }

    // If quitOnFailure is enabled and if some problems occur
    // the simulation is stopped
    if ( M_quitOnFailure && failure )
        exit( -1 );

    // If the number of iterations reaches the threshold of maxIterForReuse
    // we reset the preconditioners to force to solver to recompute it next
    // time
    if ( numIters > M_maxItersForReuse )
        resetPreconditioner();

    return numIters;
}

Real
LinearSolver::computeResidual( vector_Type& solution )
{
    if ( M_problem->getOperator() == Teuchos::null || M_problem->getRHS() == Teuchos::null )
    {
        M_displayer->leaderPrint( "SLV-  WARNING: SolverBelos can not compute the residual if the operator or the rhs is not set!\n" );
        return -1;
    }
    vector_Type res( solution.map() );
    M_problem->computeCurrResVec( &res.epetraVector(), &solution.epetraVector(), M_problem->getRHS().get() );
    Real residual;
    res.norm2( &residual );
    return residual;
}

std::string
LinearSolver::printStatus()
{
    std::ostringstream stat;
    std::string str;

    /*
     AztecOO informations:
     - Normal Convergence
     - Maximum iters reached
     - Accuracy loss
     - Ill-conditioned
     - Breakdown
     If someone has the time, he may try to find a way to report
     all these informations for the LinearSolver class as well.
     */

    if ( M_lossOfPrecision )    stat << "Accuracy loss ";
    if ( M_maxNumItersReached ) stat << "Maximum number of iterations reached ";
    str = stat.str();
    return str;
}

void
LinearSolver::setPreconditionerFromGetPot( const GetPot& dataFile, const std::string& section, PrecApplicationType precType )
{
    std::string precName = dataFile( ( section + "/prectype" ).data(), "Ifpack" );
    if ( precType == RightPreconditioner )
    {
        M_rightPreconditioner.reset( PRECFactory::instance().createObject( precName ) );
        ASSERT( M_rightPreconditioner.get() != 0, " Preconditioner not set" );
        M_rightPreconditioner->setDataFromGetPot( dataFile, section );
    }
    else
    {
        M_leftPreconditioner.reset( PRECFactory::instance().createObject( precName ) );
        ASSERT( M_leftPreconditioner.get() != 0, " Preconditioner not set" );
        M_leftPreconditioner->setDataFromGetPot( dataFile, section );
    }
}

void
LinearSolver::buildPreconditioner()
{
    LifeChrono chrono;
    Real condest( -1 );

    if ( M_leftPreconditioner )
    {
        chrono.start();
        if( !M_silent ) M_displayer->leaderPrint( "SLV-  Computing the left preconditioner...\n" );
        if ( M_baseMatrixForPreconditioner.get() == 0 )
        {
            if( !M_silent ) M_displayer->leaderPrint( "SLV-  Build preconditioner using the problem matrix\n" );
            M_leftPreconditioner->buildPreconditioner( M_matrix );
        }
        else
        {
            if( !M_silent ) M_displayer->leaderPrint( "SLV-  Build preconditioner using the base matrix provided\n" );
            M_leftPreconditioner->buildPreconditioner( M_baseMatrixForPreconditioner );
        }
        condest = M_leftPreconditioner->condest();
        Teuchos::RCP<operator_Type> leftPrec( M_leftPreconditioner->preconditioner(), false );
        Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp( new Belos::EpetraPrecOp( leftPrec ) );
        M_problem->setLeftPrec( belosPrec );
        chrono.stop();
        if( !M_silent ) M_displayer->leaderPrintMax( "SLV-  Left preconditioner computed in " , chrono.diff(), " s." );
        if( !M_silent ) M_displayer->leaderPrint( "SLV-  Estimated condition number               " , condest, "\n" );
    }
    if ( M_rightPreconditioner )
    {
        chrono.start();
        if( !M_silent ) M_displayer->leaderPrint( "SLV-  Computing the right preconditioner...\n" );
        if ( M_baseMatrixForPreconditioner.get() == 0 )
        {
            if( !M_silent ) M_displayer->leaderPrint( "SLV-  Build preconditioner using the problem matrix\n" );
            M_rightPreconditioner->buildPreconditioner( M_matrix );
        }
        else
        {
            if( !M_silent ) M_displayer->leaderPrint( "SLV-  Build preconditioner using the base matrix provided\n" );
            M_rightPreconditioner->buildPreconditioner( M_baseMatrixForPreconditioner );
        }
        condest = M_rightPreconditioner->condest();
        Teuchos::RCP<operator_Type> rightPrec( M_rightPreconditioner->preconditioner(), false );
        Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp( new Belos::EpetraPrecOp( rightPrec ) );
        M_problem->setRightPrec( belosPrec );
        chrono.stop();
        if( !M_silent ) M_displayer->leaderPrintMax( "SLV-  Right preconditioner computed in " , chrono.diff(), " s." );
        if( !M_silent ) M_displayer->leaderPrint( "SLV-  Estimated condition number               " , condest, "\n" );
    }
}

void
LinearSolver::resetPreconditioner()
{
    if ( M_leftPreconditioner )
        M_leftPreconditioner->resetPreconditioner();

    if ( M_rightPreconditioner )
        M_rightPreconditioner->resetPreconditioner();
}

bool
LinearSolver::isPreconditionerSet() const
{
    return ( M_leftPreconditioner.get() != 0 && M_leftPreconditioner->preconditionerCreated() ) ||
           ( M_rightPreconditioner.get() != 0 && M_rightPreconditioner->preconditionerCreated() ) ||
             M_problem->isLeftPrec() || M_problem->isRightPrec();
}

void
LinearSolver::showMe( std::ostream& output ) const
{
    if ( M_displayer->isLeader() )
    {
        output << "SolverBelos, parameters list:" << std::endl;
        output << "-----------------------------" << std::endl;
        output << M_parameterList << endl;
        output << "-----------------------------" << std::endl;
    }
}

// ===================================================
// Set Methods
// ===================================================
void
LinearSolver::setSolverManager( const SolverManagerType& solverManager )
{
    M_solverManagerType = solverManager;
}

void
LinearSolver::setCommunicator( const boost::shared_ptr<Epetra_Comm>& comm )
{
    M_displayer->setCommunicator( comm );
}

void LinearSolver::setMatrix( matrixPtr_Type& matrix )
{
    M_matrix = matrix;
    Teuchos::RCP<Epetra_FECrsMatrix> A = Teuchos::rcp( M_matrix->matrixPtr() );
    M_problem->setOperator( A );
}

void
LinearSolver::setOperator( Epetra_Operator& oper )
{
    M_problem->setOperator( Teuchos::rcp( &oper ) );
}

void
LinearSolver::setRightHandSide( const vector_Type& rhs )
{
    Teuchos::RCP<const vector_Type::vector_type> rhsPtr( &( rhs.epetraVector() ), false );
    M_problem->setRHS( rhsPtr );
}

void
LinearSolver::setRightHandSide( const multiVector_Type& rhs )
{
    Teuchos::RCP<const multiVector_Type> rhsPtr( &( rhs ), false );
    M_problem->setRHS( rhsPtr );
}

void
LinearSolver::setPreconditioner( preconditionerPtr_Type& preconditioner, PrecApplicationType precType )
{
    if ( precType == RightPreconditioner )
    {
        // If a right Epetra_Operator exists it must be deleted
        if ( M_problem->isRightPrec() )
        {
            M_problem->setRightPrec( Teuchos::RCP<operator_Type>( Teuchos::null ) );
        }

        M_rightPreconditioner = preconditioner;
    }
    else
    {
        // If a left Epetra_Operator exists it must be deleted
        if ( M_problem->isLeftPrec() )
        {
            M_problem->setLeftPrec( Teuchos::RCP<operator_Type>( Teuchos::null ) );
        }
        M_leftPreconditioner  = preconditioner;
    }
}

void
LinearSolver::setPreconditioner( operatorPtr_Type& preconditioner, PrecApplicationType precType )
{
    if ( precType == RightPreconditioner )
    {
        // If a right LifeV::Preconditioner exists it must be deleted
        M_rightPreconditioner.reset();

        Teuchos::RCP<operator_Type> rightPrec=Teuchos::rcp( preconditioner );
        Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp( new Belos::EpetraPrecOp( rightPrec ) );
        M_problem->setRightPrec( belosPrec );
    }
    else
    {
        // If a left LifeV::Preconditioner exists it must be deleted
        M_leftPreconditioner.reset();

        Teuchos::RCP<operator_Type> leftPrec=Teuchos::rcp( preconditioner );
        Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp( new Belos::EpetraPrecOp( leftPrec ) );
        M_problem->setLeftPrec( belosPrec );
    }
}

void
LinearSolver::setParameters( const Teuchos::ParameterList& list )
{
    M_parameterList.setParameters( list );
    M_reusePreconditioner = M_parameterList.get( "Reuse preconditioner"     , false );
    Int maxIter           = M_parameterList.get( "Maximum Iterations"       , 200 );
    M_maxItersForReuse    = M_parameterList.get( "Max Iterations For Reuse" , static_cast<Int> ( maxIter*8./10. ) );
    M_quitOnFailure       = M_parameterList.get( "Quit On Failure"          , false );
    M_silent              = M_parameterList.get( "Silent"                   , false );
}

void
LinearSolver::resetParameters()
{
    M_parameterList = Teuchos::ParameterList();
}

void
LinearSolver::setReusePreconditioner( const bool reusePreconditioner )
{
    M_reusePreconditioner = reusePreconditioner;
}

void
LinearSolver::setQuitOnFailure( const bool enable )
{
    M_quitOnFailure = enable;
}

// ===================================================
// Get Methods
// ===================================================

Int
LinearSolver::numIterations() const
{
    return M_solverManager->getNumIters();
}


Real
LinearSolver::trueResidual()
{
    if ( !M_problem->isProblemSet() )
    {
        M_displayer->leaderPrint( "SLV-  WARNING: SolverBelos can not compute the residual if the linear system is not set!\n" );
        return -1;
    }
    multiVector_Type res( *( M_problem->getRHS().get() ) );
    M_problem->computeCurrResVec( &res,M_problem->getLHS().get(), M_problem->getRHS().get() );
    Real residual;
    res.Norm2( &residual );
    return residual;
}

LinearSolver::preconditionerPtr_Type&
LinearSolver::preconditioner( PrecApplicationType precType )
{
    if ( precType == RightPreconditioner )
    {
        return M_rightPreconditioner;
    }
    return M_leftPreconditioner;
}

Teuchos::ParameterList&
LinearSolver::parametersList()
{
    return M_parameterList;
}

LinearSolver::SolverManagerPtr_Type
LinearSolver::solver()
{
    return M_solverManager;
}

boost::shared_ptr<Displayer>
LinearSolver::displayer()
{
    return M_displayer;
}

// ===================================================
// Private Methods
// ===================================================

void
LinearSolver::setupSolverManager()
{
    // If a SolverManager already exists we simply clean it!
    if ( !M_solverManager.is_null() )
    {
        M_solverManager.reset();
    }

    switch ( M_solverManagerType )
    {
        case BlockCG:
            // Create the block CG iteration
            M_solverManager = rcp( new Belos::BlockCGSolMgr<Real,multiVector_Type,operator_Type>( M_problem, rcp( &M_parameterList, false ) ) );
            break;
        case PseudoBlockCG:
            // Create the pseudo block CG iteration
            M_solverManager = rcp( new Belos::PseudoBlockCGSolMgr<Real,multiVector_Type,operator_Type>( M_problem, rcp( &M_parameterList, false) ) );
            break;
        case RCG:
            M_solverManager = rcp( new Belos::RCGSolMgr<Real,multiVector_Type,operator_Type>( M_problem, rcp( &M_parameterList, false ) ) );
            break;
        case BlockFGmres:
            M_parameterList.set( "Flexible Gmres", true );
            M_solverManager = rcp( new Belos::BlockGmresSolMgr<Real,multiVector_Type,operator_Type>( M_problem, rcp( &M_parameterList, false ) ) );
            break;
        case BlockGmres:
            M_solverManager = rcp( new Belos::BlockGmresSolMgr<Real,multiVector_Type,operator_Type>( M_problem, rcp( &M_parameterList, false ) ) );
            break;
        case PseudoBlockFGmres:
            M_parameterList.set( "Flexible Gmres", true );
            M_solverManager = rcp( new Belos::PseudoBlockGmresSolMgr<Real,multiVector_Type,operator_Type>( M_problem, rcp( &M_parameterList, false ) ) );
            break;
        case PseudoBlockGmres:
            M_solverManager = rcp( new Belos::PseudoBlockGmresSolMgr<Real,multiVector_Type,operator_Type>( M_problem, rcp( &M_parameterList, false ) ) );
            break;
        case GmresPoly:
            M_solverManager = rcp( new Belos::GmresPolySolMgr<Real,multiVector_Type,operator_Type>( M_problem, rcp( &M_parameterList,false ) ) );
            break;
        case GCRODR:
            M_solverManager = rcp( new Belos::GCRODRSolMgr<Real,multiVector_Type,operator_Type>( M_problem, rcp( &M_parameterList, false ) ) );
            break;
        case PCPG:
            M_solverManager = rcp( new Belos::PCPGSolMgr<Real,multiVector_Type,operator_Type>( M_problem, rcp( &M_parameterList, false ) ) );
            break;
        case TFQMR:
            // Create TFQMR iteration
            M_solverManager = rcp( new Belos::TFQMRSolMgr<Real,multiVector_Type,operator_Type>( M_problem, rcp( &M_parameterList, false ) ) );
            break;
    }
}

} // namespace LifeV
