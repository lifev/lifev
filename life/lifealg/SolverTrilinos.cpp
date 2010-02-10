/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
            Christoph Winkelmann <christoph.winkelmann@epfl.ch>
      Date: 2004-08-29

 Copyright (C) 2004 EPFL

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file SolverTrilinos.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-08-29
*/

#include <life/lifealg/SolverTrilinos.hpp>

namespace LifeV {

// ===================================================
// Constructors
// ===================================================
SolverTrilinos::SolverTrilinos() :
    M_prec                 (),
    M_solver               (),
    M_TrilinosParameterList(),
    M_Displayer            (),
    M_tol                  ( 0. ),
    M_maxIter              ( 0 ),
    M_maxIterForReuse      ( 0 ),
    M_reusePreconditioner                (false)
{
}

SolverTrilinos::SolverTrilinos( const Epetra_Comm& comm ) :
    M_prec                 (),
    M_solver               (),
    M_TrilinosParameterList(),
    M_Displayer            ( &comm ),
    M_tol                  ( 0. ),
    M_maxIter              ( 0 ),
    M_maxIterForReuse      ( 0 ),
    M_reusePreconditioner                (false)
{
}

// ===================================================
// Get Methods
// ===================================================
int
SolverTrilinos::NumIters() const
{
    return M_solver.NumIters();
}

int
SolverTrilinos::MaxIter() const
{
    return M_maxIter;
}


double
SolverTrilinos::TrueResidual()
{
    return M_solver.TrueResidual();
}

SolverTrilinos::prec_type&
SolverTrilinos::getPrec()
{
    return M_prec;
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

void
SolverTrilinos::getAztecStatus( double status[AZ_STATUS_SIZE] )
{
    M_solver.GetAllAztecStatus( status );
}

// ===================================================
// Set Methods
// ===================================================
void
SolverTrilinos::SetCommunicator( const Epetra_Comm& comm )
{
    M_Displayer.SetCommunicator( comm );
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
    M_maxIterForReuse = dfile( ( section + "/max_iter_reuse").data(), static_cast<int> ( M_maxIter*8./10.) );
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
    SetParameters( false );
}

void
SolverTrilinos::SetParameters(bool cerr_warning_if_unused)
{
    M_solver.SetParameters(M_TrilinosParameterList,cerr_warning_if_unused);
}

void
SolverTrilinos::setTolMaxiter(const double tol, const int maxiter)
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


// ===================================================
// Methods
// ===================================================
int
SolverTrilinos::solve( vector_type& x, const vector_type& b )
{
    M_solver.SetLHS( &x.getEpetraVector() );
    // The Solver from Aztecoo takes a non const (because of rescaling?)
    // We should be careful if you use scaling
    Epetra_FEVector* bVectorPtr ( const_cast<Epetra_FEVector*> (&b.getEpetraVector()) );
    M_solver.SetRHS( bVectorPtr );

    int    maxiter(M_maxIter);
    double mytol  (M_tol);
    int status;

    if ( isPrecSet() && M_prec->precType().compare("AztecOO") )
        M_solver.SetPrecOperator(M_prec->getPrec());

    status = M_solver.Iterate(maxiter, mytol);

#ifdef DEBUG

     M_Displayer.comm().Barrier();
     M_Displayer.leaderPrint( "  o-  Number of iterations = ", M_solver.NumIters());
     M_Displayer.leaderPrint( "  o-  Norm of the true residual = ", M_solver.TrueResidual());
     M_Displayer.leaderPrint( "  o-  Norm of the true ratio    = ",  M_solver.ScaledResidual());
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
         int olditer = M_solver.NumIters();
         status = M_solver.Iterate(maxiter, mytol);

#ifdef DEBUG
         M_Displayer.comm().Barrier();
         M_Displayer.leaderPrint( "  o-  Second run: number of iterations = ", M_solver.NumIters());
         M_Displayer.leaderPrint( "  o-  Norm of the true residual = ",  M_solver.TrueResidual());
         M_Displayer.leaderPrint( "  o-  Norm of the true ratio    = ",  M_solver.ScaledResidual());
#endif
         return(M_solver.NumIters() + olditer);
     }

     return(M_solver.NumIters());
}

double
SolverTrilinos::computeResidual( vector_type& x, vector_type& b )
{
    vector_type Ax(x.getMap());
    vector_type res(b);

    M_solver.GetUserMatrix()->Apply(x.getEpetraVector(),Ax.getEpetraVector());

    res.getEpetraVector().Update(1, Ax.getEpetraVector(), -1);

    double residual;

    res.Norm2(&residual);

    return residual;
}

std::string
SolverTrilinos::printStatus()
{

    std::ostringstream stat;
    std::string str;

    double status[AZ_STATUS_SIZE];
    getAztecStatus( status );

    if( status[AZ_why] == AZ_normal         ) stat << "Normal Convergence    ";
    else if( status[AZ_why] == AZ_maxits    ) stat << "Maximum iters reached ";
    else if( status[AZ_why] == AZ_loss      ) stat << "Accuracy loss         ";
    else if( status[AZ_why] == AZ_ill_cond  ) stat << "Ill-conditioned       ";
    else if( status[AZ_why] == AZ_breakdown ) stat << "Breakdown             ";

    stat << setw(12) << "res = " << status[AZ_scaled_r];
    stat << setw(4)  << " " << (int)status[AZ_its] << " iters. ";
    stat << std::endl;

    str = stat.str();
    return str;
}

bool
SolverTrilinos::isPrecSet() const
{
    return ( M_prec.get() !=0 && M_prec->getPrec() != 0 );
}

void SolverTrilinos::buildPreconditioner( matrix_ptrtype& prec)
{
    Chrono chrono;
    double condest(-1);

    chrono.start();

    M_Displayer.leaderPrint("      Computing the precond ...                ");

    M_prec->buildPreconditioner(prec);

    condest = M_prec->Condest();
    chrono.stop();

    M_Displayer.leaderPrintMax( "done in " , chrono.diff() );
    M_Displayer.leaderPrint("      Estimated condition number = " , condest );
}

void SolverTrilinos::setUpPrec(const GetPot& dataFile,  const std::string& section)
{
    std::string precType = dataFile( (section + "/prectype").data(), "Ifpack");
    M_prec.reset( PRECFactory::instance().createObject( precType ) );

    ASSERT(M_prec.get() != 0, " Preconditioner not set");
    M_prec->setSolver( *this );
    M_prec->setDataFromGetPot( dataFile, section );
}


int SolverTrilinos::solveSystem( const vector_type&      rhsFull,
                                 vector_type&      sol,
                                 matrix_ptrtype&   baseMatrixForPreconditioner)

{

    bool retry(true);

    Chrono chrono;

    M_Displayer.leaderPrint("      Setting up the solver ...                \n");

    if (baseMatrixForPreconditioner.get() == 0)
        M_Displayer.leaderPrint("      Warning: baseMatrixForPreconditioner is empty     \n");

    if ( !isPrecSet() || !M_reusePreconditioner  )
    {
        buildPreconditioner(baseMatrixForPreconditioner);
        // do not retry if I am recomputing the preconditioner
        retry = false;
    }
    else
    {
        M_Displayer.leaderPrint("      Reusing  precond ...                 \n");
    }

    int numIter = solveSystem(rhsFull, sol, M_prec);

    // If we do not want to retry, return now.
    // otherwise rebuild the preconditioner and solve again:
    if ( numIter < 0  && retry )
    {
        chrono.start();

        M_Displayer.leaderPrint("     Iterative solver failed, numiter = " , - numIter);
        M_Displayer.leaderPrint("     maxIterSolver = " , M_maxIter );
        M_Displayer.leaderPrint("     retrying:          ");

        buildPreconditioner(baseMatrixForPreconditioner);

        chrono.stop();
        M_Displayer.leaderPrintMax( "done in " , chrono.diff() );
        // Solving again, but only once (retry = false)
        numIter = solveSystem(rhsFull, sol, M_prec);

        if (numIter < 0)
            M_Displayer.leaderPrint(" ERROR: Iterative solver failed again.\n");
    }

    if ( abs(numIter) > M_maxIterForReuse)
        precReset();

    return numIter;
}

int SolverTrilinos::solveSystem( const vector_type&  rhsFull,
                                 vector_type&        sol,
                                 prec_type&          prec )

{
    M_Displayer.leaderPrint("      Solving system ...                   \n");

    setPreconditioner(prec);

    Chrono chrono;
    chrono.start();
    int numIter = solve( sol, rhsFull );
    chrono.stop();
    M_Displayer.leaderPrintMax( "       ... done in " , chrono.diff() );

    if (numIter >= M_maxIter)
        numIter = -numIter;

    return numIter;
}

} // namespace LifeV
