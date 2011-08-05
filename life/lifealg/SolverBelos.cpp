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

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include "Teuchos_RCPBoostSharedPtrConversions.hpp"

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/LifeV.hpp>
#include <life/lifealg/SolverBelos.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
SolverBelos::SolverBelos() :
        M_leftPreconditioner   (),
        M_rightPreconditioner  (),
        M_solverManager        (),
        M_problem              ( new LinearProblem_type() ),
        M_parameterList        (),
        M_displayer            ( new Displayer() ),
        M_maxIterForReuse      ( 0 ),
        M_reusePreconditioner  (false)
{
}

SolverBelos::SolverBelos( const boost::shared_ptr<Epetra_Comm>& comm ) :
        M_leftPreconditioner   (),
        M_rightPreconditioner  (),
        M_solverManager        (),
        M_problem              ( new LinearProblem_type() ),
        M_parameterList        (),
        M_displayer            ( new Displayer(comm) ),
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
    M_problem->setLHS( rcp(&solution.epetraVector()) );
    // WARNING: This remark is stated in Aztec00. I do not know if
    //          it is also valid for Belos
    //        > The Solver from Aztecoo takes a non const (because of rescaling?)
    //        > We should be careful if you use scaling
    Epetra_FEVector* rhsVectorPtr ( const_cast<Epetra_FEVector*> (&rhs.epetraVector()) );
    M_problem->setRHS( rcp(rhsVectorPtr) );

    //Int  maxiter(M_maxIter);
    //Real mytol  (M_tolerance);
    Int status;

    //if ( isPreconditionerSet() && M_rightPreconditioner->preconditionerType().compare("AztecOO") )
        //M_solver.SetPrecOperator(M_rightPreconditioner->preconditioner()); //Left prec??????????????????????????????? todo

    //status = M_solver.Iterate(maxiter, mytol);

#ifdef HAVE_LIFEV_DEBUG
    M_displayer->comm()->Barrier();
    //M_displayer->leaderPrint( "  o-  Number of iterations = ", M_solver.NumIters());
    //M_displayer->leaderPrint( "  o-  Norm of the true residual = ", M_solver.TrueResidual());
    //M_displayer->leaderPrint( "  o-  Norm of the true ratio    = ",  M_solver.ScaledResidual());
#endif

    /* try to solve again (reason may be:
      -2 "Aztec status AZ_breakdown: numerical breakdown"
      -3 "Aztec status AZ_loss: loss of precision"
      -4 "Aztec status AZ_ill_cond: GMRES hessenberg ill-conditioned"
    */
    if ( status <= -2 )
    {
        //maxiter     = M_maxIter; //todo
        //mytol       = M_tolerance; //todo
        //Int oldIter = M_solver.NumIters();
        //status      = M_solver.Iterate(maxiter, mytol);

#ifdef HAVE_LIFEV_DEBUG
        M_displayer->comm()->Barrier();
        //M_displayer->leaderPrint( "  o-  Second run: number of iterations = ", M_solver.NumIters());
        //M_displayer->leaderPrint( "  o-  Norm of the true residual = ",  M_solver.TrueResidual());
        //M_displayer->leaderPrint( "  o-  Norm of the true ratio    = ",  M_solver.ScaledResidual());
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

    M_problem->getOperator()->Apply( solution.epetraVector(), Ax.epetraVector() );

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
    // todo deal with preconditioner
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

    Int numIter = solveSystem( rhsFull, solution, M_rightPreconditioner );

    // If we do not want to retry, return now.
    // otherwise rebuild the preconditioner and solve again:
    if ( numIter < 0  && retry )
    {
        chrono.start();

        M_displayer->leaderPrint( "SLV-  Iterative solver failed, numiter = " , - numIter );
        //M_displayer->leaderPrint( "SLV-  maxIterSolver = " , M_maxIter );
        M_displayer->leaderPrint( "SLV-  retrying:          " );

        buildPreconditioner( baseMatrixForPreconditioner );

        chrono.stop();
        M_displayer->leaderPrintMax( "done in " , chrono.diff() );
        // Solving again, but only once (retry = false)
        numIter = solveSystem( rhsFull, solution, M_rightPreconditioner );

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
    M_rightPreconditioner.reset( PRECFactory::instance().createObject( precType ) );

    ASSERT( M_rightPreconditioner.get() != 0, " Preconditioner not set" );
    //M_rightPreconditioner->setSolver( *this ); GWENOL TODO
    M_rightPreconditioner->setDataFromGetPot( dataFile, section );
}

void SolverBelos::buildPreconditioner( matrix_ptrtype& preconditioner )
{
    LifeChrono chrono;
    Real condest(-1);

    chrono.start();

    M_displayer->leaderPrint( "SLV-  Computing the precond ...                " );

    M_rightPreconditioner->buildPreconditioner( preconditioner ); //todo

    condest = M_rightPreconditioner->condest();
    chrono.stop();

    M_displayer->leaderPrintMax( "done in " , chrono.diff() );
    M_displayer->leaderPrint( "SLV-  Estimated condition number               " , condest, "\n" );
}

void SolverBelos::resetPreconditioner()
{
    M_rightPreconditioner->resetPreconditioner();
    M_leftPreconditioner->resetPreconditioner();
}

bool
SolverBelos::isPreconditionerSet() const
{
    return ( M_leftPreconditioner.get() !=0 && M_leftPreconditioner->preconditionerCreated() )||( M_rightPreconditioner.get() !=0 && M_rightPreconditioner->preconditionerCreated() );
}

void
SolverBelos::showMe( std::ostream& output ) const
{
    output << "SolverBelos, parameters list:" << std::endl;
    output << "-----------------------------" << std::endl;
    output << M_parameterList << endl;
    output << "-----------------------------" << std::endl;
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
    /*
    M_matrix = matrix.matrixPtr();
    boost::shared_ptr<const Epetra_FECrsMatrix> const_matrix=boost::const_pointer_cast<const Epetra_FECrsMatrix>(M_matrix);
    boost::shared_ptr<const Epetra_Operator> const_operator=boost::static_pointer_cast<const Epetra_Operator>(const_matrix);
    std::cout << "pointer: " << M_matrix.get() << std::endl;
    //M_sA = Teuchos::rcp(M_matrix);
    M_sA = Teuchos::rcp(const_operator);
    std::cout << "pointer: " << M_sA.get() << std::endl;

    M_problem->setOperator( M_sA );
    */

    M_matrix = matrix.matrixPtr();
    //std::cout << "pointer: " << M_matrix.get() << std::endl;
    Teuchos::RCP<Epetra_FECrsMatrix> M_sA = Teuchos::rcp(M_matrix);
    //std::cout << "pointer: " << M_sA.get() << std::endl;

    M_problem->setOperator( M_sA );
}

void
SolverBelos::setOperator( Epetra_Operator& oper )
{
    M_problem->setOperator(rcp(&oper));
}

void
SolverBelos::setPreconditioner( prec_type& preconditioner, PrecApplicationType precType )
{
    if(precType == RightPreconditioner)
    {
        M_rightPreconditioner = preconditioner;
    }
    else
    {
        M_leftPreconditioner  = preconditioner;
    }
}

void
SolverBelos::setPreconditioner( comp_prec_type& preconditioner, PrecApplicationType precType )
{
    //M_solver.SetPrecOperator( preconditioner.get() );
}

void
SolverBelos::setParameters( const GetPot& dataFile, const std::string& section )
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
    Real tolerance = dataFile( ( section + "/tol" ).data(), 1.e-6 );
    M_parameterList.set( "Convergence Tolerance", tolerance );

    // Maximum Number of iterations
    Int maxIter         = dataFile( ( section + "/max_iter"      ).data(), 200 );
    M_parameterList.set( "Maximum Iterations", maxIter );

    M_maxIterForReuse = dataFile( ( section + "/max_iter_reuse").data(), static_cast<Int> ( maxIter*8./10.) );
    M_reusePreconditioner = dataFile( (section + "/reuse").data(), M_reusePreconditioner );



    // GMRES PARAMETERS

    // Krylov space dimension
    //M_parameterList.set( "kspace", dataFile( ( section + "/kspace" ).data(), M_maxIter ) );

    // Gram-Schmidt algorithm
    //M_parameterList.set( "orthog", dataFile( ( section + "/orthog" ).data(), AZ_classic ) );

    // r-vector
    //M_parameterList.set( "aux_vec", dataFile( ( section + "/aux_vec" ).data(), AZ_resid ) );
}

void
SolverBelos::setParameters( const Teuchos::ParameterList& list )
{
    M_parameterList.setParameters( list );
}

void
SolverBelos::resetParameters()
{
    M_parameterList = Teuchos::ParameterList();
}

void
SolverBelos::setReusePreconditioner( const bool reusePreconditioner )
{
    M_reusePreconditioner = reusePreconditioner;
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


Real
SolverBelos::trueResidual()
{
    //return M_solver.TrueResidual();
    return 1.0;
}

SolverBelos::prec_type&
SolverBelos::preconditioner( PrecApplicationType precType )
{
    if(precType == RightPreconditioner)
    {
        return M_rightPreconditioner;
    }
    return M_leftPreconditioner;
}

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

boost::shared_ptr<Displayer>
SolverBelos::displayer()
{
    return M_displayer;
}

void
SolverBelos::setupSolverManager()
{
    // Create the block CG iteration
    //M_solverManager = rcp( new Belos::BlockCGSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false)) );

    // Create the block GMRes iteration
    // Create the flexible, block GMRes iteration
    M_solverManager = rcp( new Belos::BlockGmresSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false) ) );

    /*
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
    */
}

} // namespace LifeV
