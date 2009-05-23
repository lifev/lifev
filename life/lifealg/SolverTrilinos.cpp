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

#include <life/lifecore/debug.hpp>

#include <life/lifecore/GetPot.hpp>

#include <life/lifealg/SolverTrilinos.hpp>
#include "Epetra_Comm.h"

#include <iomanip>

namespace LifeV
{
// namespace Epetra
// {


SolverTrilinos::SolverTrilinos(Epetra_Comm& comm)
    :
    M_prec  (),
    M_solver(),
    M_TrilinosParameterList(),
    M_Displayer(comm)
{
}
SolverTrilinos::SolverTrilinos()
    :
    M_prec  (),
    M_solver(),
    M_TrilinosParameterList(),
    M_Displayer()
{
}


int
SolverTrilinos::NumIters()
{
    return M_solver.NumIters();
}


double
SolverTrilinos::TrueResidual()
{
    return M_solver.TrueResidual();
}

void SolverTrilinos::setMatrix( matrix_type& m)
{
    M_solver.SetUserMatrix(&m.getEpetraMatrix());
}

void SolverTrilinos::setOperator( Epetra_Operator& op)
{
    M_solver.SetUserOperator(&op);
}

//! set Epetra_Operator preconditioner
void SolverTrilinos::setPreconditioner( prec_type _prec )
{
    M_prec = _prec;
}


void SolverTrilinos::SetParameters(bool cerr_warning_if_unused)
{
    M_solver.SetParameters(M_TrilinosParameterList,cerr_warning_if_unused);
}


void SolverTrilinos::useGMRES(const int restart)
{
    M_TrilinosParameterList.set("solver", "AZ_gmres");
    //M_TrilinosParameterList.set("solver", "AZ_gmres_condnum");
    M_TrilinosParameterList.set("kspace", restart);
}

void SolverTrilinos::useCG()
{
    M_TrilinosParameterList.set("solver", "AZ_cg");
    //M_TrilinosParameterList.set("solver", "AZ_cg_condnum");
}


void SolverTrilinos::useBICGSTAB()
{
    M_TrilinosParameterList.set("solver", "AZ_bicgstab");
    //M_TrilinosParameterList.set("solver", "AZ_cg_condnum");
}

void SolverTrilinos::useJacobyPrec()
{
    M_TrilinosParameterList.set("precond", "Jacobi");
}

void SolverTrilinos::useNoPrec()
{
    M_TrilinosParameterList.set("precond", "none");
}

void SolverTrilinos::useDDPrec()
{
    M_TrilinosParameterList.set("precond", "dom_decomp");
    //M_TrilinosParameterList.set("subdomain_solve", AZ_ilu );
    //M_TrilinosParameterList.set("graph_fill", 2 );
    M_TrilinosParameterList.set("drop", 0 );
    M_TrilinosParameterList.set("ilut_fill", 10 );
    M_TrilinosParameterList.set("overlap", 1 );
    M_TrilinosParameterList.set("type_overlap", "symmetric" );
}

void SolverTrilinos::useNeumannPrec()
{
    M_TrilinosParameterList.set("precond", "Neumann");
    M_TrilinosParameterList.set("poly_ord", 10 );
}

void SolverTrilinos::setReuse()
{
    M_TrilinosParameterList.set("pre_calc", "reuse");
}

void SolverTrilinos::setDataFromGetPot( const GetPot& dfile, const std::string& section )
{
    //
    // solver
    //
    M_TrilinosParameterList.set("solver",
                             dfile( ( section + "/solver" ).data(), "gmres" ));
    //
    // scaling
    //
    M_TrilinosParameterList.set("scaling",
                             dfile( ( section + "/scaling" ).data(), "none" ));
    //
    // precond
    //
    M_TrilinosParameterList.set("precond",
                             dfile( ( section + "/precond" ).data(), "none" ));
    //
    // Convergence type
    //
    M_TrilinosParameterList.set("conv",
                             dfile( ( section + "/conv" ).data(), "rhs" ));
    //
    // output
    //
    M_TrilinosParameterList.set("output",
                             dfile( ( section + "/output" ).data(), "all_res" ));

    //
    // other variables

    M_TrilinosParameterList.set("pre_calc",
                             dfile( ( section + "/pre_calc" ).data(), AZ_calc ));


    M_maxIter = dfile( ( section + "/max_iter" ).data(), 200 );
    M_TrilinosParameterList.set("max_iter", M_maxIter);

    M_TrilinosParameterList.set("poly_ord",
                             dfile( ( section + "/poly_ord" ).data(), 3 ));

    M_TrilinosParameterList.set("overlap",
                             dfile( ( section + "/overlap" ).data(), AZ_none ));

    M_TrilinosParameterList.set("type_overlap",
                             dfile( ( section + "/type_overlap" ).data(), AZ_standard ));

    M_TrilinosParameterList.set("kspace",
                             dfile( ( section + "/kspace" ).data(), M_maxIter ));

    M_TrilinosParameterList.set("orthog",
                             dfile( ( section + "/orthog" ).data(), AZ_classic ));

    M_TrilinosParameterList.set("aux_vec",
                             dfile( ( section + "/aux_vec" ).data(), AZ_resid ));

    M_TrilinosParameterList.set("reorder",
                             dfile( ( section + "/reorder" ).data(), 1 ));

    M_TrilinosParameterList.set("keep_info",
                             dfile( ( section + "/keep_info" ).data(), 0 ));

    //
    // subdomain solver
    //

    M_TrilinosParameterList.set("subdomain_solve",
                             dfile( ( section + "/subdomain_solve" ).data(),
                                    "ilut" ));

    M_TrilinosParameterList.set("graph_fill",
                             dfile( ( section + "/graph_fill" ).data(), 0 ));

    M_TrilinosParameterList.set("init_guess",
                             dfile( ( section + "/init_guess" ).data(), AZ_NOT_ZERO ));

    M_TrilinosParameterList.set("keep_kvecs",
                             dfile( ( section + "/keep_kvecs" ).data(), 0 ));

    M_TrilinosParameterList.set("apply_kvecs",
                             dfile( ( section + "/apply_kvecs" ).data(), AZ_FALSE ));

    M_TrilinosParameterList.set("orth_kvecs",
                             dfile( ( section + "/orth_kvecs" ).data(), AZ_FALSE ));

    M_TrilinosParameterList.set("ignore_scaling",
                             dfile( ( section + "/ignore_scaling" ).data(), AZ_FALSE ));

    M_TrilinosParameterList.set("check_update_size",
                             dfile( ( section + "/check_update_size" ).data(),
                                    AZ_FALSE ));
    //


    M_tol = dfile( ( section + "/tol" ).data(), 1.0e-06 );
    M_TrilinosParameterList.set("tol", M_tol);


//     M_TrilinosParameterList.set("drop",
//                              dfile( ( section + "/drop" ).data(), 0.0 ));

//     M_TrilinosParameterList.set("ilut_fill",
//                              dfile( ( section + "/ilut_fill" ).data(), 1. ));

    M_TrilinosParameterList.set("omega",
                             dfile( ( section + "/omega" ).data(), 1. ));

    M_TrilinosParameterList.set("update_reduction",
                             dfile( ( section + "/update_reduction" ).data(), 10e10 ));

    M_maxIterSolver   = dfile(( section + "/max_iter").data(), -1);

}


void SolverTrilinos::setTolMaxiter(const double tol, const int maxiter)
{
    if (tol > 0)
        {
            M_tol = tol;
            M_TrilinosParameterList.set("tol", M_tol);
        }

    if (maxiter >= 0)
        {
            M_maxIter = maxiter;
            M_TrilinosParameterList.set("max_iter", M_maxIter);
        }

}

void  SolverTrilinos::SetVerbose(const VerboseLevel verb)
{
    switch (verb) {
    case NONE:
        M_TrilinosParameterList.set("output", "none");
        return;
    case SUMMARY:
        M_TrilinosParameterList.set("output", "summary");
        return;
    case LAST:
    default:
        M_TrilinosParameterList.set("output", "last");
        return;
    }
}


int
SolverTrilinos::solve( vector_type& x, vector_type& b )
{

    M_solver.SetLHS(&x.getEpetraVector());
    M_solver.SetRHS(&b.getEpetraVector());

    SetParameters(true);

//    M_TrilinosParameterList.set("kspace", 100);
//    M_TrilinosParameterList.TrilinosParameterListShowMe();

    int    maxiter(M_maxIter);
    double mytol  (M_tol);
    int status;
     if ( precSet() )
//         {
        M_solver.SetPrecOperator(M_prec->getPrec());

     status = M_solver.Iterate(maxiter, mytol);
     //        }
//     else
//         {
//              status = M_solver.AdaptiveIterate(maxiter, 3, mytol);
//         }
    /* if status:
       0  AZ_normal
       1  AZ_maxits
       if < 0 see AztecOO.cpp
    */

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


int SolverTrilinos::solveSystem(  vector_type&      rhsFull,
                                  vector_type&      sol,
                                  matrix_ptrtype&   basePrecMatrix,
                                  bool const        reuse,
                                  bool const        retry)
{
    Chrono chrono;
    M_Displayer.leaderPrint("      Setting up the solver ...                \n");

    double condest(-1);
    if ( !M_prec->set() || !reuse  )
    {
        chrono.start();

        M_Displayer.leaderPrint("      Computing the precond ...                \n");

        M_prec->buildPreconditioner(basePrecMatrix);

        condest = M_prec->Condest();

        setPreconditioner(M_prec);

        chrono.stop();
        M_Displayer.leaderPrintMax( "done in " , chrono.diff() );
        M_Displayer.leaderPrint("      Estimated condition number = " , condest );
    }
    else
    {
        M_Displayer.leaderPrint("      Reusing  precond ...                \n");
    }

    M_Displayer.leaderPrint("      Solving system ...                       \n");

    chrono.start();
    int numIter = solve(sol, rhsFull);
    chrono.stop();
    M_Displayer.leaderPrintMax( "       ... done in " , chrono.diff() );

    // If we do not want to retry, return now.
    // otherwise rebuild the preconditioner and solve again:
    if (numIter >= M_maxIterSolver && retry)
    {
        chrono.start();

        M_Displayer.leaderPrint("     Iterative solver failed, numiter = " , numIter);
        M_Displayer.leaderPrint("     maxIterSolver = " , M_maxIterSolver );
        M_Displayer.leaderPrint("     retrying: rebuilding prec ...          ");

        if (basePrecMatrix.get())
        {
            M_prec->precReset();

            M_prec->buildPreconditioner(basePrecMatrix);

            condest = M_prec->Condest();

            setPreconditioner(M_prec);
        }
        chrono.stop();
        M_Displayer.leaderPrintMax( "done in " , chrono.diff() );
        M_Displayer.leaderPrint("     Estimated condition number = " , condest );
        // Solving again, but only once (retry = false)
        numIter = solveSystem(rhsFull, sol,
                              basePrecMatrix, reuse, false);

        if (numIter >= M_maxIterSolver)
            M_Displayer.leaderPrint(" ERROR: Iterative solver failed again.\n");
        return -numIter;
    }
    return numIter;
}

void SolverTrilinos::setUpPrec(const GetPot& dataFile,  const std::string& section)
{
    std::string precType = dataFile( (section + "/prectype").data(), "Ifpack");
    M_prec.reset( PRECFactory::instance().createObject( precType ) );
    ASSERT(M_prec.get() != 0, "Oseen : Preconditioner not set");
    M_prec->setDataFromGetPot( dataFile, section );
}

void SolverTrilinos::setPrec(prec_raw_type* prec)
{M_prec.reset(prec);}

} // namespace LifeV

