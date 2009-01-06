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

#include "life/lifealg/SolverTrilinos.hpp"
#include "Epetra_Comm.h"



namespace LifeV
{
// namespace Epetra
// {


SolverTrilinos::SolverTrilinos()
    :
    M_prec  (),
    M_solver(),
    M_TrilinosParameterList()
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
                             dfile( ( section + "/kspace" ).data(), 100 ));

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

    if ( precSet() )
        M_solver.SetPrecOperator(M_prec->getPrec());

    int    maxiter(M_maxIter);
    double mytol  (M_tol);

    int status = M_solver.Iterate(maxiter, mytol);
    /* if status:
       0  AZ_normal
       1  AZ_maxits
       if < 0 see AztecOO.cpp
    */

//#ifdef DEBUG
    const Epetra_Comm* Comm(0);
    if (M_solver.GetUserOperator()) Comm = &M_solver.GetUserOperator()->Comm();
    else if (M_solver.GetUserMatrix()) Comm = &M_solver.GetUserMatrix()->Comm();

    if (! Comm)
        {
            Comm->Barrier();

            if( Comm->MyPID() == 0 ) {
                std::cout << "  o-  Solver performed " << M_solver.NumIters()
                          << " iterations.\n";
                std::cout << "  o-  Norm of the true residual = " << M_solver.TrueResidual() << std::endl;
                std::cout << "  o-  Norm of the true ratio    = " << M_solver.ScaledResidual() << std::endl;
            }
        }
//#endif

    if (status == -2 || status == -3 ) // try to solve again (reason may be:
        {
            maxiter = M_maxIter;
            mytol = M_tol;
            int olditer = M_solver.NumIters();
            status = M_solver.Iterate(maxiter, mytol);

//#ifdef DEBUG
    if (! Comm)
        {
            Comm->Barrier();

            if( Comm->MyPID() == 0 ) {
                std::cout << "  o-  Second run: solver performed " << M_solver.NumIters()
                          << " iterations.\n";
                std::cout << "  o-  Norm of the true residual = " << M_solver.TrueResidual() << std::endl;
                std::cout << "  o-  Norm of the true ratio    = " << M_solver.ScaledResidual() << std::endl;
            }
        }
//#endif
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



// } //namespace Epetra




} // namespace LifeV

