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
    @brief SolverTrilinos

    @author Simone Deparis   <simone.deparis@epfl.ch>
    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 08-11-2006
 */

#include <life/lifecore/life.hpp>
#include <life/lifealg/SolverTrilinos.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
SolverTrilinos::SolverTrilinos() :
        M_preconditioner       (),
        M_solver               (),
        M_TrilinosParameterList(),
        M_displayer            ( new Displayer() ),
        M_tolerance            ( 0. ),
        M_maxIter              ( 0 ),
        M_maxIterForReuse      ( 0 ),
        M_reusePreconditioner  (false)
{
}

SolverTrilinos::SolverTrilinos( const boost::shared_ptr<Epetra_Comm>& comm ) :
        M_preconditioner       (),
        M_solver               (),
        M_TrilinosParameterList(),
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
SolverTrilinos::solve( vector_type& solution, const vector_type& rhs )
{
    M_solver.SetLHS( &solution.getEpetraVector() );
    // The Solver from Aztecoo takes a non const (because of rescaling?)
    // We should be careful if you use scaling
    Epetra_FEVector* rhsVectorPtr ( const_cast<Epetra_FEVector*> (&rhs.getEpetraVector()) );
    M_solver.SetRHS( rhsVectorPtr );

    Int  maxiter(M_maxIter);
    Real mytol  (M_tolerance);
    Int status;

    if ( isPrecSet() && M_preconditioner->precType().compare("AztecOO") )
        M_solver.SetPrecOperator(M_preconditioner->getPrec());

    status = M_solver.Iterate(maxiter, mytol);

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
        Int oldIter = M_solver.NumIters();
        status      = M_solver.Iterate(maxiter, mytol);

#ifdef HAVE_LIFEV_DEBUG
        M_displayer->comm()->Barrier();
        M_displayer->leaderPrint( "  o-  Second run: number of iterations = ", M_solver.NumIters());
        M_displayer->leaderPrint( "  o-  Norm of the true residual = ",  M_solver.TrueResidual());
        M_displayer->leaderPrint( "  o-  Norm of the true ratio    = ",  M_solver.ScaledResidual());
#endif
        return( M_solver.NumIters() + oldIter );
    }

    return( M_solver.NumIters() );
}

Real
SolverTrilinos::computeResidual( vector_type& solution, vector_type& rhs )
{
    vector_type Ax ( solution.getMap() );
    vector_type res( rhs );

    M_solver.GetUserMatrix()->Apply( solution.getEpetraVector(), Ax.getEpetraVector() );

    res.getEpetraVector().Update( 1, Ax.getEpetraVector(), -1 );

    Real residual;

    res.Norm2( &residual );

    return residual;
}

std::string
SolverTrilinos::printStatus()
{

    std::ostringstream stat;
    std::string str;

    Real status[AZ_STATUS_SIZE];
    getAztecStatus( status );

    if ( status[AZ_why] == AZ_normal         ) stat << "Normal Convergence    ";
    else if ( status[AZ_why] == AZ_maxits    ) stat << "Maximum iters reached ";
    else if ( status[AZ_why] == AZ_loss      ) stat << "Accuracy loss         ";
    else if ( status[AZ_why] == AZ_ill_cond  ) stat << "Ill-conditioned       ";
    else if ( status[AZ_why] == AZ_breakdown ) stat << "Breakdown             ";

    stat << setw(12) << "res = " << status[AZ_scaled_r];
    stat << setw(4)  << " " << (Int)status[AZ_its] << " iters. ";
    stat << std::endl;

    str = stat.str();
    return str;
}

Int SolverTrilinos::solveSystem( const vector_type& rhsFull,
                                 vector_type&       solution,
                                 matrix_ptrtype&    baseMatrixForPreconditioner )

{

    bool retry( true );

    Chrono chrono;

    M_displayer->leaderPrint( "      Setting up the solver ...                \n" );

    if ( baseMatrixForPreconditioner.get() == 0 )
        M_displayer->leaderPrint( "      Warning: baseMatrixForPreconditioner is empty     \n" );

    if ( !isPrecSet() || !M_reusePreconditioner  )
    {
        buildPreconditioner( baseMatrixForPreconditioner );
        // do not retry if I am recomputing the preconditioner
        retry = false;
    }
    else
    {
        M_displayer->leaderPrint( "      Reusing precond ...                 \n" );
    }

    Int numIter = solveSystem( rhsFull, solution, M_preconditioner );

    // If we do not want to retry, return now.
    // otherwise rebuild the preconditioner and solve again:
    if ( numIter < 0  && retry )
    {
        chrono.start();

        M_displayer->leaderPrint( "     Iterative solver failed, numiter = " , - numIter );
        M_displayer->leaderPrint( "     maxIterSolver = " , M_maxIter );
        M_displayer->leaderPrint( "     retrying:          " );

        buildPreconditioner( baseMatrixForPreconditioner );

        chrono.stop();
        M_displayer->leaderPrintMax( "done in " , chrono.diff() );
        // Solving again, but only once (retry = false)
        numIter = solveSystem( rhsFull, solution, M_preconditioner );

        if ( numIter < 0 )
            M_displayer->leaderPrint( " ERROR: Iterative solver failed again.\n" );
    }

    if ( abs(numIter) > M_maxIterForReuse )
        precReset();

    return numIter;
}

void SolverTrilinos::setUpPrec( const GetPot& dataFile,  const std::string& section )
{
    std::string precType = dataFile( (section + "/prectype").data(), "Ifpack" );
    M_preconditioner.reset( PRECFactory::instance().createObject( precType ) );

    ASSERT( M_preconditioner.get() != 0, " Preconditioner not set" );
    M_preconditioner->setSolver( *this );
    M_preconditioner->setDataFromGetPot( dataFile, section );
}

void SolverTrilinos::buildPreconditioner( matrix_ptrtype& preconditioner )
{
    Chrono chrono;
    Real condest(-1);

    chrono.start();

    M_displayer->leaderPrint( "      Computing the precond ...                " );

    M_preconditioner->buildPreconditioner( preconditioner );

    condest = M_preconditioner->Condest();
    chrono.stop();

    M_displayer->leaderPrintMax( "done in " , chrono.diff() );
    M_displayer->leaderPrint( "      Estimated condition number               " , condest, "\n" );
}

void SolverTrilinos::precReset()
{
    M_preconditioner->precReset();
}

bool
SolverTrilinos::isPrecSet() const
{
    return ( M_preconditioner.get() !=0 && M_preconditioner->preconditionerCreated() );
}

void
SolverTrilinos::showMe( std::ostream& output ) const
{
    output << "showMe must be implemented for the SolverTrilinos class" << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void
SolverTrilinos::setCommunicator( const boost::shared_ptr<Epetra_Comm>& comm )
{
    M_displayer->setCommunicator( comm );
}

void SolverTrilinos::setMatrix( matrix_type& matrix )
{
    M_matrix = matrix.getMatrixPtr();
    M_solver.SetUserMatrix( M_matrix.get() );
}

void
SolverTrilinos::setOperator( Epetra_Operator& oper )
{
    M_solver.SetUserOperator( &oper );
}

void
SolverTrilinos::setPreconditioner( prec_type& preconditioner )
{
    M_preconditioner = preconditioner;
}

void
SolverTrilinos::setPreconditioner( comp_prec_type& preconditioner )
{
    M_solver.SetPrecOperator( preconditioner.get() );
}

void
SolverTrilinos::setDataFromGetPot( const GetPot& dataFile, const std::string& section )
{
    // SOLVER PARAMETERS

    // Solver type
    M_TrilinosParameterList.set( "solver",  dataFile( ( section + "/solver" ).data(), "gmres" ) );

    // Residual expression
    M_TrilinosParameterList.set( "conv",    dataFile( ( section + "/conv" ).data(), "rhs" ) );

    // Scaling
    M_TrilinosParameterList.set( "scaling", dataFile( ( section + "/scaling" ).data(), "none" ) );

    // Output
    M_TrilinosParameterList.set( "output",  dataFile( ( section + "/output" ).data(), "all" ) );

    // Tolerance
    M_tolerance = dataFile( ( section + "/tol" ).data(), 1.e-6 );
    M_TrilinosParameterList.set( "tol", M_tolerance );

    // Maximum Number of iterations
    M_maxIter         = dataFile( ( section + "/max_iter"      ).data(), 200 );
    M_maxIterForReuse = dataFile( ( section + "/max_iter_reuse").data(), static_cast<Int> ( M_maxIter*8./10.) );
    M_reusePreconditioner = dataFile( (section + "/reuse").data(), M_reusePreconditioner );

    M_TrilinosParameterList.set( "max_iter", M_maxIter );

    // GMRES PARAMETERS

    // Krylov space dimension
    M_TrilinosParameterList.set( "kspace", dataFile( ( section + "/kspace" ).data(), M_maxIter ) );

    // Gram-Schmidt algorithm
    M_TrilinosParameterList.set( "orthog", dataFile( ( section + "/orthog" ).data(), AZ_classic ) );

    // r-vector
    M_TrilinosParameterList.set( "aux_vec", dataFile( ( section + "/aux_vec" ).data(), AZ_resid ) );


    // SET PARAMETERS
    setParameters( false );
}

void
SolverTrilinos::setParameters( bool cerrWarningIfUnused )
{
    M_solver.SetParameters( M_TrilinosParameterList, cerrWarningIfUnused );
}

void
SolverTrilinos::setTolMaxiter( const Real tolerance, const Int maxIter )
{
    if ( tolerance > 0 )
    {
        M_tolerance = tolerance;
        M_TrilinosParameterList.set( "tol", M_tolerance );
    }

    if ( maxIter >= 0 )
    {
        M_maxIter = maxIter;
        M_TrilinosParameterList.set( "max_iter", M_maxIter );
    }
}

void
SolverTrilinos::setReusePreconditioner( const bool reusePreconditioner )
{
    M_reusePreconditioner = reusePreconditioner;
}

boost::shared_ptr<Displayer>
SolverTrilinos::displayer()
{
    return M_displayer;
}

// ===================================================
// Get Methods
// ===================================================
Int
SolverTrilinos::NumIters() const
{
    return M_solver.NumIters();
}

Int
SolverTrilinos::MaxIter() const
{
    return M_maxIter;
}


Real
SolverTrilinos::TrueResidual()
{
    return M_solver.TrueResidual();
}

SolverTrilinos::prec_type&
SolverTrilinos::getPrec()
{
    return M_preconditioner;
}

void
SolverTrilinos::getAztecStatus( Real status[AZ_STATUS_SIZE] )
{
    M_solver.GetAllAztecStatus( status );
}

Teuchos::ParameterList&
SolverTrilinos::getParameterList()
{
    return M_TrilinosParameterList;
}

AztecOO&
SolverTrilinos::getSolver()
{
    return M_solver;
}

} // namespace LifeV
