/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file SolverAmesos.cpp

*/




/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Gilles Fourestey <christophe.prudhomme@epfl.ch>
      Date: 2009-06-09

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
   \file SolverAmesos.cpp
   \author Gilles Fourestey <gilles.fourestey@epfl.ch>
   \date 2004-08-29
*/

#include <life/lifecore/debug.hpp>

#include <life/lifecore/GetPot.hpp>

#include <life/lifealg/SolverAmesos.hpp>
#include "Epetra_Comm.h"

#include <iomanip>

namespace LifeV
{
// namespace Epetra
// {


SolverAmesos::SolverAmesos(Epetra_Comm& comm)
    :
    M_matrix               (),
    M_factory              (),
    M_problem              (),
    M_TrilinosParameterList(),
    M_redistribute         (true),
    M_printTiming          (false),
    M_printStatus          (false),
    M_Displayer            (&comm),
    M_comm                 (comm)
{
}

// SolverAmesos::SolverAmesos():
//     M_prec                 (),
//     M_factory              (),
//     M_TrilinosParameterList(),
//     M_Displayer            (),

// {
//     assert(false);
// }


int
SolverAmesos::NumIters()
{
    return 1;
}


double
SolverAmesos::TrueResidual()
{
    return 0.;//M_factory.TrueResidual();
}

void SolverAmesos::setMatrix(matrix_type& m)
{
    M_matrix = m.getMatrixPtr();
    M_problem.SetOperator(M_matrix.get());
}

void SolverAmesos::setOperator( Epetra_Operator& op)
{
    ASSERT(false,"SolverAmesos::setOperator: not coded");
}


void SolverAmesos::SetParameters(bool cerr_warning_if_unused)
{

}



void SolverAmesos::setDataFromGetPot( const GetPot& dfile, const std::string& section )
{
    M_solverType    = dfile(( section + "/amesos/solvertype").data(), "Superludist");
    M_redistribute  = dfile(( section + "/amesos/redistribute").data(), true);
    M_printTiming   = dfile(( section + "/amesos/print_timing").data(), false);
    M_printStatus   = dfile(( section + "/amesos/print_status").data(), false);
}


void SolverAmesos::setTolMaxiter(const double tol, const int maxiter)
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

void  SolverAmesos::SetVerbose(const VerboseLevel verb)
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


// int
// SolverAmesos::solve( vector_type& x, vector_type& b )
// {



//     //M_factory.SetLHS(&x.getEpetraVector());
//     //M_factory.SetRHS(&b.getEpetraVector());
//     //M_Displayer.comm().Barrier();

//     //SetParameters(true);

// //    M_TrilinosParameterList.set("kspace", 100);
// //    M_TrilinosParameterList.TrilinosParameterListShowMe();

//     //int    maxiter(M_maxIter);
//     //double mytol  (M_tol);
//     //int status;


//     if ( precSet() )
// //         {
//          M_factory.SetPrecOperator(M_prec->getPrec());

//      status = M_factory.Iterate(maxiter, mytol);
//      //        }
// //     else
// //         {
// //              status = M_factory.AdaptiveIterate(maxiter, 3, mytol);
// //         }
//     /* if status:
//        0  AZ_normal
//        1  AZ_maxits
//        if < 0 see AztecOO.cpp
//     */

// #ifdef DEBUG

//      M_Displayer.comm().Barrier();
//      M_Displayer.leaderPrint( "  o-  Number of iterations = ", M_factory.NumIters());
//      M_Displayer.leaderPrint( "  o-  Norm of the true residual = ", M_factory.TrueResidual());
//      M_Displayer.leaderPrint( "  o-  Norm of the true ratio    = ",  M_factory.ScaledResidual());
// #endif

//      /* try to solve again (reason may be:
//        -2 "Aztec status AZ_breakdown: numerical breakdown"
//        -3 "Aztec status AZ_loss: loss of precision"
//        -4 "Aztec status AZ_ill_cond: GMRES hessenberg ill-conditioned"
//      */
//      if (status <= -2 )
//      {
//          maxiter = M_maxIter;
//          mytol = M_tol;
//          int olditer = M_factory.NumIters();
//          status = M_factory.Iterate(maxiter, mytol);

// #ifdef DEBUG
//          M_Displayer.comm().Barrier();
//          M_Displayer.leaderPrint( "  o-  Second run: number of iterations = ", M_factory.NumIters());
//          M_Displayer.leaderPrint( "  o-  Norm of the true residual = ",  M_factory.TrueResidual());
//          M_Displayer.leaderPrint( "  o-  Norm of the true ratio    = ",  M_factory.ScaledResidual());
// #endif
//          return(M_factory.NumIters() + olditer);
//      }

//      return(M_factory.NumIters());

// }

double
SolverAmesos::computeResidual( vector_type& x, vector_type& b )
{
    vector_type Ax(x.getMap());
    vector_type res(b);

    res.getEpetraVector().Update(1, Ax.getEpetraVector(), -1);

    double residual;

    res.Norm2(&residual);

    return residual;
}



std::string
SolverAmesos::printStatus()
{


    // 1) The symbolic factorization
    //    (parameter doesn't always exist)
    std::cout << "  Amesos: Total symbolic factorization time " << M_sfact_time << std::endl;

    // 2) The numeric factorization
    //    (always exists if NumericFactorization() is called)
    std::cout << "  Amesos: Total numeric factorization time  " << M_nfact_time << std::endl;
    // 3) Solving the linear system
    //    (always exists if Solve() is called)
    std::cout << "  Amesos: Total solve time                  " << M_solve_time << std::endl;

    // 4) Converting the matrix to the accepted format for the solver
    //    (always exists if SymbolicFactorization() is called)
    std::cout << "  Amesos: matrix convertion time            " << M_mtx_conv_time << std::endl;

    // 5) Redistributing the matrix for each solve to the accepted format for the solver
    std::cout << "  Amesos: Total matrix redistribution time  " << M_mtx_redist_time << std::endl;

    // 6) Redistributing the vector for each solve to the accepted format for the solver
    std::cout << "  Amesos: Total vector redistribution time  " << M_vec_redist_time << std::endl;
}


int SolverAmesos::solveSystem(  vector_type&      rhsFull,
                                vector_type&      sol,
                                matrix_ptrtype&   basePrecMatrix,
                                bool const        reuse,
                                bool const        retry)
{

    if (M_comm.MyPID() == 0)
         std::cout << "       Solving the system ... " << std::endl;

    Amesos_BaseSolver* Solver;

    M_problem.SetLHS(&sol.getEpetraVector());
    M_problem.SetRHS(&rhsFull.getEpetraVector());

    Solver = M_factory.Create(M_solverType, M_problem);

    if (Solver == 0) {
        std::cerr << std::endl  << std::endl;
        std::cerr << "SolverAmesos: Selected solver << " << M_solverType << " is not available. Bailing out." << std::endl;
        // return ok not to break the test harness
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        exit(EXIT_SUCCESS);
    }

    Teuchos::ParameterList List;

    List.set("PrintTiming",  M_printTiming);
    List.set("PrintStatus",  M_printStatus);
    List.set("Redistribute", M_redistribute);

    Solver->SetParameters(List);

    if (M_comm.MyPID() == 0)
        std::cout << "          Starting symbolic factorization  ..." << std::flush;
    Solver->SymbolicFactorization();
    if (M_comm.MyPID() == 0)
        std::cout << " ok." << std::endl;

    //std::cout << "m shape = " << Solver->MatrixShapeOK() << std::endl;
  // you can change the matrix values here
    if (M_comm.MyPID() == 0)
        std::cout << "          Starting numeric factorization   ..." << std::flush;
    Solver->NumericFactorization();
    if (M_comm.MyPID() == 0)
        std::cout << " ok." << std::endl;

    // you can change LHS and RHS here
    if (M_comm.MyPID() == 0)
        std::cout << "          Starting solution phase          ..." << std::flush;
    Solver->Solve();
    if (M_comm.MyPID() == 0)
        std::cout << " ok." << std::endl;

    // you can get the timings here
    Teuchos::ParameterList TimingsList;
    Solver->GetTiming( TimingsList );

    if (M_printTiming) Solver->PrintTiming();
    if (M_printStatus) Solver->PrintStatus();
    //printStatus();

    return 0;
}

void SolverAmesos::setUpPrec(const GetPot& dataFile,  const std::string& section)
{
    return;
}
//     //     std::string precType = dataFile( (section + "/prectype").data(), "Ifpack");

//     //     M_prec.reset( PRECFactory::instance().createObject( precType ) );
//     //     ASSERT(M_prec.get() != 0, "Oseen : Preconditioner not set");
//     //     M_prec->setDataFromGetPot( dataFile, section );
// }

// void SolverAmesos::setPrec(prec_raw_type* prec)
// {
//     return;
// }

} // namespace LifeV

