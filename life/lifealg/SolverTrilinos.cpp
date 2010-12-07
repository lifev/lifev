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

#include <lifeconfig.h>
#include <life/lifealg/SolverTrilinos.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
SolverTrilinos::SolverTrilinos() :
        M_prec                 (),
        M_solver               (),
        M_TrilinosParameterList(),
        M_displayer            (new Displayer()),
        M_tol                  ( 0. ),
        M_maxIter              ( 0 ),
        M_maxIterForReuse      ( 0 ),
        M_reusePreconditioner                (false)
{
}

SolverTrilinos::SolverTrilinos( const boost::shared_ptr<Epetra_Comm>& comm ) :
        M_prec                 (),
        M_solver               (),
        M_TrilinosParameterList(),
        M_displayer            ( new Displayer(comm) ),
        M_tol                  ( 0. ),
        M_maxIter              ( 0 ),
        M_maxIterForReuse      ( 0 ),
        M_reusePreconditioner                (false)
{
}

// ===================================================
// Methods
// ===================================================
Int
SolverTrilinos::solve( vector_type& x, const vector_type& b )
{
    M_solver.SetLHS( &x.getEpetraVector() );
    // The Solver from Aztecoo takes a non const (because of rescaling?)
    // We should be careful if you use scaling
    Epetra_FEVector* bVectorPtr ( const_cast<Epetra_FEVector*> (&b.getEpetraVector()) );
    M_solver.SetRHS( bVectorPtr );

    Int    maxiter(M_maxIter);
    Real mytol  (M_tol);
    Int status;

    if ( isPrecSet() && M_prec->precType().compare("AztecOO") )
        M_solver.SetPrecOperator(M_prec->getPrec());

    status = M_solver.Iterate(maxiter, mytol);

#ifdef DEBUG

    M_displayer->comm().Barrier();
    M_displayer->leaderPrint( "  o-  Number of iterations = ", M_solver.NumIters());
    M_displayer->leaderPrint( "  o-  Norm of the true residual = ", M_solver.TrueResidual());
    M_displayer->leaderPrint( "  o-  Norm of the true ratio    = ",  M_solver.ScaledResidual());
#endif

    /* try to solve again (reason may be:
      -2 "Aztec status AZ_breakdown: numerical breakdown"
      -3 "Aztec status AZ_loss: loss of precision"
      -4 "Aztec status AZ_ill_cond: GMRES hessenberg ill-conditioned"
    */
    if (status <= -2 )
    {
        maxiter = M_maxIter;
        mytol = M_tol;
        Int olditer = M_solver.NumIters();
        status = M_solver.Iterate(maxiter, mytol);

#ifdef DEBUG
        M_displayer->comm().Barrier();
        M_displayer->leaderPrint( "  o-  Second run: number of iterations = ", M_solver.NumIters());
        M_displayer->leaderPrint( "  o-  Norm of the true residual = ",  M_solver.TrueResidual());
        M_displayer->leaderPrint( "  o-  Norm of the true ratio    = ",  M_solver.ScaledResidual());
#endif
        return(M_solver.NumIters() + olditer);
    }

    return(M_solver.NumIters());
}

Real
SolverTrilinos::computeResidual( vector_type& x, vector_type& b )
{
    vector_type Ax(x.getMap());
    vector_type res(b);

    M_solver.GetUserMatrix()->Apply(x.getEpetraVector(),Ax.getEpetraVector());

    res.getEpetraVector().Update(1, Ax.getEpetraVector(), -1);

    Real residual;

    res.Norm2(&residual);

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

Int SolverTrilinos::solveSystem( const vector_type&      rhsFull,
                                 vector_type&      sol,
                                 matrix_ptrtype&   baseMatrixForPreconditioner)

{

    bool retry(true);

    Chrono chrono;

    M_displayer->leaderPrint("      Setting up the solver ...                \n");

    if (baseMatrixForPreconditioner.get() == 0)
        M_displayer->leaderPrint("      Warning: baseMatrixForPreconditioner is empty     \n");

    if ( !isPrecSet() || !M_reusePreconditioner  )
    {
        buildPreconditioner(baseMatrixForPreconditioner);
        // do not retry if I am recomputing the preconditioner
        retry = false;
    }
    else
    {
        M_displayer->leaderPrint("      Reusing precond ...                 \n");
    }

    Int numIter = solveSystem(rhsFull, sol, M_prec);

    // If we do not want to retry, return now.
    // otherwise rebuild the preconditioner and solve again:
    if ( numIter < 0  && retry )
    {
        chrono.start();

        M_displayer->leaderPrint("     Iterative solver failed, numiter = " , - numIter);
        M_displayer->leaderPrint("     maxIterSolver = " , M_maxIter );
        M_displayer->leaderPrint("     retrying:          ");

        buildPreconditioner(baseMatrixForPreconditioner);

        chrono.stop();
        M_displayer->leaderPrintMax( "done in " , chrono.diff() );
        // Solving again, but only once (retry = false)
        numIter = solveSystem(rhsFull, sol, M_prec);

        if (numIter < 0)
            M_displayer->leaderPrint(" ERROR: Iterative solver failed again.\n");
    }

    if ( abs(numIter) > M_maxIterForReuse)
        precReset();

    return numIter;
}

void SolverTrilinos::setUpPrec(const GetPot& dataFile,  const std::string& section)
{
    std::string precType = dataFile( (section + "/prectype").data(), "Ifpack");
    M_prec.reset( PRECFactory::instance().createObject( precType ) );

    ASSERT(M_prec.get() != 0, " Preconditioner not set");
    M_prec->setSolver( *this );
    M_prec->setDataFromGetPot( dataFile, section );
}

void SolverTrilinos::buildPreconditioner( matrix_ptrtype& prec)
{
    Chrono chrono;
    Real condest(-1);

    chrono.start();

    M_displayer->leaderPrint("      Computing the precond ...                ");

    M_prec->buildPreconditioner(prec);

    condest = M_prec->Condest();
    chrono.stop();

    M_displayer->leaderPrintMax( "done in " , chrono.diff() );
    M_displayer->leaderPrint("      Estimated condition number               " , condest, "\n" );
}

void SolverTrilinos::precReset()
{
    M_prec->precReset();
}

bool
SolverTrilinos::isPrecSet() const
{
    return ( M_prec.get() !=0 && M_prec->preconditionerCreated() );
}

// ===================================================
// Set Methods
// ===================================================
void
SolverTrilinos::setCommunicator( const boost::shared_ptr<Epetra_Comm>& comm )
{
    M_displayer->setCommunicator( comm );
}

void SolverTrilinos::setMatrix( matrix_type& m)
{
    M_matrix = m.getMatrixPtr();
    M_solver.SetUserMatrix(M_matrix.get());
}

void
SolverTrilinos::setOperator( Epetra_Operator& op)
{
    M_solver.SetUserOperator( &op );
}

void
SolverTrilinos::setPreconditioner( prec_type& _prec )
{
    M_prec = _prec;
}

void
SolverTrilinos::setPreconditioner( comp_prec_type& _prec )
{
    M_solver.SetPrecOperator(_prec.get());
}

void
SolverTrilinos::setDataFromGetPot( const GetPot& dfile, const std::string& section )
{
    // SOLVER PARAMETERS

    // Solver type
    M_TrilinosParameterList.set("solver",  dfile( ( section + "/solver" ).data(), "gmres" ));

    // Residual expression
    M_TrilinosParameterList.set("conv",    dfile( ( section + "/conv" ).data(), "rhs" ));

    // Scaling
    M_TrilinosParameterList.set("scaling", dfile( ( section + "/scaling" ).data(), "none" ));

    // Output
    M_TrilinosParameterList.set("output",  dfile( ( section + "/output" ).data(), "all" ));

    // Tolerance
    M_tol = dfile( ( section + "/tol" ).data(), 1.e-6 );
    M_TrilinosParameterList.set("tol", M_tol);

    // Maximum Number of iterations
    M_maxIter         = dfile( ( section + "/max_iter"      ).data(), 200 );
    M_maxIterForReuse = dfile( ( section + "/max_iter_reuse").data(), static_cast<Int> ( M_maxIter*8./10.) );
    M_reusePreconditioner = dfile( (section + "/reuse").data(), M_reusePreconditioner);

    M_TrilinosParameterList.set("max_iter", M_maxIter);

    // GMRES PARAMETERS

    // Krylov space dimension
    M_TrilinosParameterList.set("kspace", dfile( ( section + "/kspace" ).data(), M_maxIter ));

    // Gram-Schmidt algorithm
    M_TrilinosParameterList.set("orthog", dfile( ( section + "/orthog" ).data(), AZ_classic ));

    // r-vector
    M_TrilinosParameterList.set("aux_vec", dfile( ( section + "/aux_vec" ).data(), AZ_resid ));


    // SET PARAMETERS
    setParameters( false );
}

void
SolverTrilinos::setParameters(bool cerr_warning_if_unused)
{
    M_solver.SetParameters(M_TrilinosParameterList,cerr_warning_if_unused);
}

void
SolverTrilinos::setTolMaxiter(const Real tol, const Int maxiter)
{
    if ( tol > 0 )
    {
        M_tol = tol;
        M_TrilinosParameterList.set("tol", M_tol);
    }

    if ( maxiter >= 0 )
    {
        M_maxIter = maxiter;
        M_TrilinosParameterList.set("max_iter", M_maxIter);
    }
}

void
SolverTrilinos::setReusePreconditioner( const bool reuse )
{
    M_reusePreconditioner = reuse;
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
    return M_prec;
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
