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
    @brief SolverBelos

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 03-08-2011
 */

#include <life/lifecore/LifeV.hpp>
#include <life/lifealg/SolverBelos.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
SolverBelos::SolverBelos() :
        M_preconditioner       (),
        M_solverManager        (),
        M_problem              (),
        M_parameterList        (),
        M_displayer            ( new Displayer() ),
        M_tolerance            ( 0. ),
        M_maxIter              ( 0 ),
        M_maxIterForReuse      ( 0 ),
        M_reusePreconditioner  (false)
{
}

SolverBelos::SolverBelos( const boost::shared_ptr<Epetra_Comm>& comm ) :
        M_preconditioner       (),
        M_solverManager        (),
        M_problem              (),
        M_parameterList        (),
        M_displayer            ( new Displayer(comm) ),
        M_tolerance            ( 0. ),
        M_maxIter              ( 0 ),
        M_maxIterForReuse      ( 0 ),
        M_reusePreconditioner  (false)
{
}

// ===================================================
// Methods
// ===================================================
Int
SolverBelos::solve( vector_type& solution, const vector_type& rhs )
{
    //M_solver.SetLHS( &solution.epetraVector() );
    // The Solver from Aztecoo takes a non const (because of rescaling?)
    // We should be careful if you use scaling
    Epetra_FEVector* rhsVectorPtr ( const_cast<Epetra_FEVector*> (&rhs.epetraVector()) );
    //M_solver.SetRHS( rhsVectorPtr );

    Int  maxiter(M_maxIter);
    Real mytol  (M_tolerance);
    Int status;

    //if ( isPreconditionerSet() && M_preconditioner->preconditionerType().compare("AztecOO") )
        //M_solver.SetPrecOperator(M_preconditioner->preconditioner());

    //status = M_solver.Iterate(maxiter, mytol);

#ifdef HAVE_LIFEV_DEBUG
    M_displayer->comm()->Barrier();
    M_displayer->leaderPrint( "  o-  Number of iterations = ", M_solver.NumIters());
    M_displayer->leaderPrint( "  o-  Norm of the true residual = ", M_solver.TrueResidual());
    M_displayer->leaderPrint( "  o-  Norm of the true ratio    = ",  M_solver.ScaledResidual());
#endif

    /* try to solve again (reason may be:
      -2 "Aztec status AZ_breakdown: numerical breakdown"
      -3 "Aztec status AZ_loss: loss of precision"
      -4 "Aztec status AZ_ill_cond: GMRES hessenberg ill-conditioned"
    */
    if ( status <= -2 )
    {
        maxiter     = M_maxIter;
        mytol       = M_tolerance;
        //Int oldIter = M_solver.NumIters();
        //status      = M_solver.Iterate(maxiter, mytol);

#ifdef HAVE_LIFEV_DEBUG
        M_displayer->comm()->Barrier();
        M_displayer->leaderPrint( "  o-  Second run: number of iterations = ", M_solver.NumIters());
        M_displayer->leaderPrint( "  o-  Norm of the true residual = ",  M_solver.TrueResidual());
        M_displayer->leaderPrint( "  o-  Norm of the true ratio    = ",  M_solver.ScaledResidual());
#endif
        //return( M_solver.NumIters() + oldIter );
    }

    //return( M_solver.NumIters() );
    return 0; // GWENOL todo
}

Real
SolverBelos::computeResidual( vector_type& solution, vector_type& rhs )
{
    vector_type Ax ( solution.map() );
    vector_type res( rhs );

    //M_solver.GetUserMatrix()->Apply( solution.epetraVector(), Ax.epetraVector() );

    res.epetraVector().Update( 1, Ax.epetraVector(), -1 );

    Real residual;

    res.norm2( &residual );

    return residual;
}

std::string
SolverBelos::printStatus()
{

    std::ostringstream stat;
    std::string str;

    /*
    Real status[AZ_STATUS_SIZE];
    aztecStatus( status );

    if ( status[AZ_why] == AZ_normal         ) stat << "Normal Convergence    ";
    else if ( status[AZ_why] == AZ_maxits    ) stat << "Maximum iters reached ";
    else if ( status[AZ_why] == AZ_loss      ) stat << "Accuracy loss         ";
    else if ( status[AZ_why] == AZ_ill_cond  ) stat << "Ill-conditioned       ";
    else if ( status[AZ_why] == AZ_breakdown ) stat << "Breakdown             ";

    stat << setw(12) << "res = " << status[AZ_scaled_r];
    stat << setw(4)  << " " << (Int)status[AZ_its] << " iters. ";
    stat << std::endl;
    */

    str = stat.str();
    return str;
}

Int SolverBelos::solveSystem( const vector_type& rhsFull,
                                 vector_type&       solution,
                                 matrix_ptrtype&    baseMatrixForPreconditioner )

{

    bool retry( true );

    LifeChrono chrono;

    M_displayer->leaderPrint( "SLV-  Setting up the solver ...                \n" );

    if ( baseMatrixForPreconditioner.get() == 0 )
        M_displayer->leaderPrint( "SLV-  Warning: baseMatrixForPreconditioner is empty     \n" );

    if ( !isPreconditionerSet() || !M_reusePreconditioner  )
    {
        buildPreconditioner( baseMatrixForPreconditioner );
        // do not retry if I am recomputing the preconditioner
        retry = false;
    }
    else
    {
        M_displayer->leaderPrint( "SLV-  Reusing precond ...                 \n" );
    }

    Int numIter = solveSystem( rhsFull, solution, M_preconditioner );

    // If we do not want to retry, return now.
    // otherwise rebuild the preconditioner and solve again:
    if ( numIter < 0  && retry )
    {
        chrono.start();

        M_displayer->leaderPrint( "SLV-  Iterative solver failed, numiter = " , - numIter );
        M_displayer->leaderPrint( "SLV-  maxIterSolver = " , M_maxIter );
        M_displayer->leaderPrint( "SLV-  retrying:          " );

        buildPreconditioner( baseMatrixForPreconditioner );

        chrono.stop();
        M_displayer->leaderPrintMax( "done in " , chrono.diff() );
        // Solving again, but only once (retry = false)
        numIter = solveSystem( rhsFull, solution, M_preconditioner );

        if ( numIter < 0 )
            M_displayer->leaderPrint( " ERROR: Iterative solver failed again.\n" );
    }

    if ( std::abs(numIter) > M_maxIterForReuse )
        resetPreconditioner();

    return numIter;
}

void SolverBelos::setupPreconditioner( const GetPot& dataFile,  const std::string& section )
{
    std::string precType = dataFile( (section + "/prectype").data(), "Ifpack" );
    M_preconditioner.reset( PRECFactory::instance().createObject( precType ) );

    ASSERT( M_preconditioner.get() != 0, " Preconditioner not set" );
    //M_preconditioner->setSolver( *this ); GWENOL TODO
    M_preconditioner->setDataFromGetPot( dataFile, section );
}

void SolverBelos::buildPreconditioner( matrix_ptrtype& preconditioner )
{
    LifeChrono chrono;
    Real condest(-1);

    chrono.start();

    M_displayer->leaderPrint( "SLV-  Computing the precond ...                " );

    M_preconditioner->buildPreconditioner( preconditioner );

    condest = M_preconditioner->condest();
    chrono.stop();

    M_displayer->leaderPrintMax( "done in " , chrono.diff() );
    M_displayer->leaderPrint( "SLV-  Estimated condition number               " , condest, "\n" );
}

void SolverBelos::resetPreconditioner()
{
    M_preconditioner->resetPreconditioner();
}

bool
SolverBelos::isPreconditionerSet() const
{
    return ( M_preconditioner.get() !=0 && M_preconditioner->preconditionerCreated() );
}

void
SolverBelos::showMe( std::ostream& output ) const
{
    output << "showMe must be implemented for the SolverBelos class" << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void
SolverBelos::setCommunicator( const boost::shared_ptr<Epetra_Comm>& comm )
{
    M_displayer->setCommunicator( comm );
}

void SolverBelos::setMatrix( matrix_type& matrix )
{
    M_matrix = matrix.matrixPtr();
    //M_solver.SetUserMatrix( M_matrix.get() );
}

void
SolverBelos::setOperator( Epetra_Operator& oper )
{
    //M_solver.SetUserOperator( &oper );
}

void
SolverBelos::setPreconditioner( prec_type& preconditioner )
{
    M_preconditioner = preconditioner;
}

void
SolverBelos::setPreconditioner( comp_prec_type& preconditioner )
{
    //M_solver.SetPrecOperator( preconditioner.get() );
}

void
SolverBelos::setDataFromGetPot( const GetPot& dataFile, const std::string& section )
{
    // SOLVER PARAMETERS

    // Solver type
    //M_parameterList.set( "solver",  dataFile( ( section + "/solver" ).data(), "gmres" ) );

    // Residual expression
    //M_parameterList.set( "conv",    dataFile( ( section + "/conv" ).data(), "rhs" ) );

    // Scaling
    //M_parameterList.set( "scaling", dataFile( ( section + "/scaling" ).data(), "none" ) );

    // Output
    //M_parameterList.set( "output",  dataFile( ( section + "/output" ).data(), "all" ) );

    // Tolerance
    M_tolerance = dataFile( ( section + "/tol" ).data(), 1.e-6 );
    M_parameterList.set( "Convergence Tolerance", M_tolerance );

    // Maximum Number of iterations
    M_maxIter         = dataFile( ( section + "/max_iter"      ).data(), 200 );
    M_maxIterForReuse = dataFile( ( section + "/max_iter_reuse").data(), static_cast<Int> ( M_maxIter*8./10.) );
    M_reusePreconditioner = dataFile( (section + "/reuse").data(), M_reusePreconditioner );

    M_parameterList.set( "Maximum Iterations", M_maxIter );

    // GMRES PARAMETERS

    // Krylov space dimension
    //M_parameterList.set( "kspace", dataFile( ( section + "/kspace" ).data(), M_maxIter ) );

    // Gram-Schmidt algorithm
    //M_parameterList.set( "orthog", dataFile( ( section + "/orthog" ).data(), AZ_classic ) );

    // r-vector
    //M_parameterList.set( "aux_vec", dataFile( ( section + "/aux_vec" ).data(), AZ_resid ) );


    // SET PARAMETERS
    setParameters( false );
}

void
SolverBelos::setParameters( bool cerrWarningIfUnused )
{
    //M_solver.SetParameters( M_parameterList, cerrWarningIfUnused );
}

void
SolverBelos::setTolerance( const Real tolerance )
{
    if ( tolerance > 0 )
    {
        M_tolerance = tolerance;
        M_parameterList.set( "Convergence Tolerance", M_tolerance );
    }
}

void
SolverBelos::setMaxNumIterations( const Int maxIter )
{
    if ( maxIter >= 0 )
    {
        M_maxIter = maxIter;
        M_parameterList.set( "Maximum Iterations", M_maxIter );
    }
}

void
SolverBelos::setReusePreconditioner( const bool reusePreconditioner )
{
    M_reusePreconditioner = reusePreconditioner;
}

boost::shared_ptr<Displayer>
SolverBelos::displayer()
{
    return M_displayer;
}

// ===================================================
// Get Methods
// ===================================================
//!
Int
SolverBelos::numIterations() const
{
    return M_solverManager->getNumIters();
}

Int
SolverBelos::maxNumIterations() const
{
    return M_maxIter;
}


Real
SolverBelos::trueResidual()
{
    //return M_solver.TrueResidual();
    return 1.0;
}

SolverBelos::prec_type&
SolverBelos::preconditioner()
{
    return M_preconditioner;
}

/*
void
SolverBelos::aztecStatus( Real status[AZ_STATUS_SIZE] )
{
    //M_solver.GetAllAztecStatus( status );
}
*/

Teuchos::ParameterList&
SolverBelos::getParametersList()
{
    return M_parameterList;
}

SolverBelos::SolverManager_ptrtype
SolverBelos::solver()
{
    return M_solverManager;
}

void
SolverBelos::setupSolverManager()
{
    // Create the block CG iteration
    M_solverManager = rcp( new Belos::BlockCGSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false)) );

    // Create the block GMRes iteration
    // Create the flexible, block GMRes iteration
    M_solverManager = rcp( new Belos::BlockGmresSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false) ) );

    M_solverManager = rcp( new Belos::GCRODRSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false) ) );

    M_solverManager = rcp( new Belos::GmresPolySolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false) ) );

    M_solverManager = rcp( new Belos::PCPGSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false)) );

    // Create the pseudo block CG iteration
    M_solverManager = rcp( new Belos::PseudoBlockCGSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false)) );

    // Create the pseudo block GMRes iteration
    M_solverManager = rcp( new Belos::PseudoBlockGmresSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false) ) );

    M_solverManager = rcp( new Belos::RCGSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false)) );

    // Create TFQMR iteration
    M_solverManager = rcp( new Belos::TFQMRSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false)) );
}

void
SolverBelos::createEpetraProblem()
{
    //M_problem = rcp( new LinearProblem_type( A, LHS, RHS ) );
}

} // namespace LifeV
