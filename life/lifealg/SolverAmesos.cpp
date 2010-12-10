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
    @brief Solver Amesos

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 29-08-2004
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_Comm.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifeconfig.h>
#include <life/lifealg/SolverAmesos.hpp>
#include <life/lifecore/debug.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifecore/GetPot.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
SolverAmesos::SolverAmesos( const comm_PtrType& comm ) :
        M_matrix               (),
        M_problem              (),
        M_solver               (),
        M_trilinosParameterList(),
        M_displayer            ( comm )
{
}

// ===================================================
// Methods
// ===================================================
Real
SolverAmesos::computeResidual( const vector_type& solution, const vector_type& rhs )
{
    vector_type Ax(solution.getMap());
    vector_type res(rhs);

    res.getEpetraVector().Update(1, Ax.getEpetraVector(), -1);

    Real residual;

    res.Norm2(&residual);

    return residual;
}

Int
SolverAmesos::solveSystem( vector_type&    rhsFull,
                           vector_type&    solution,
                           matrix_ptrtype& /*basePrecMatrix*/ )
{
    M_displayer.leaderPrint( "      Amesos solving system ...                " );
    Chrono chrono;
    chrono.start();

    M_problem.SetLHS( &solution.getEpetraVector() );
    M_problem.SetRHS( &rhsFull.getEpetraVector() );

    M_solver->Solve();

    chrono.stop();
    M_displayer.leaderPrintMax( "done in " , chrono.diff() );

    return 0;
}

void
SolverAmesos::printStatus()
{
    /*
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
    */

    if ( M_trilinosParameterList.get( "PrintTiming", false ) )
        M_solver->PrintTiming();

    if ( M_trilinosParameterList.get( "PrintStatus", false ) )
        M_solver->PrintStatus();
}

bool SolverAmesos::isPrecSet() const
{
    return true;
}

void SolverAmesos::precReset()
{

}

void SolverAmesos::setUpPrec( const GetPot& /*dataFile*/, const std::string& /*section*/ )
{

}

void SolverAmesos::setReusePreconditioner( const bool& /*reusePreconditioner*/ )
{

}

void SolverAmesos::showMe( std::ostream& output ) const
{
    output << "showMe must be implemented for the SolverAmesos class" << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void SolverAmesos::setMatrix( const matrix_type& matrix )
{
    M_matrix = matrix.getMatrixPtr();
    M_problem.SetOperator( M_matrix.get() );

    // After setting the matrix we can perform symbolic & numeric factorization
    M_solver->SymbolicFactorization();
    M_solver->NumericFactorization();
}

void SolverAmesos::setOperator( const Epetra_Operator& /*oper*/ )
{
    ASSERT( false, "SolverAmesos::setOperator: not coded" );
}

void SolverAmesos::setDataFromGetPot( const GetPot& dataFile, const std::string& section )
{
    // Status parameters
    M_trilinosParameterList.set( "OutputLevel",  dataFile( ( section + "/amesos/outputlevel").data(), 0 ) );
    M_trilinosParameterList.set( "PrintStatus",  dataFile( ( section + "/amesos/print_status").data(), false ) );
    M_trilinosParameterList.set( "PrintTiming",  dataFile( ( section + "/amesos/print_timing").data(), false ) );
    M_trilinosParameterList.set( "ComputeVectorNorms", dataFile( ( section + "/amesos/computevectornorms").data(), false ) );
    M_trilinosParameterList.set( "ComputeTrueResidual", dataFile( ( section + "/amesos/computeresidual").data(), false ) );

    // Control parameters
    M_trilinosParameterList.set( "AddZeroToDiag",  dataFile( ( section + "/amesos/addzerotodiag").data(), false ) );
    M_trilinosParameterList.set( "Refactorize", dataFile( ( section + "/amesos/refactorize").data(), false ) );
    M_trilinosParameterList.set( "RcondThreshold", dataFile( ( section + "/amesos/rcondthreshold").data(), 1.e-2) );
    M_trilinosParameterList.set( "Redistribute", dataFile( ( section + "/amesos/redistribute").data(), true ) ); // SuperLU
    M_trilinosParameterList.set( "MaxProcs", dataFile( ( section + "/amesos/maxprocs").data(), -1) ); // ScalaPack
    M_trilinosParameterList.set( "ScaleMethod", dataFile( ( section + "/amesos/scalemethod").data(), 1) );

    // Type of the matrix: symmetric, SDP, general
    M_trilinosParameterList.set( "MatrixProperty", dataFile( ( section + "/amesos/matrixproperty").data(), "general" ) );

    // Create the solver
    createSolver( dataFile( ( section + "/amesos/solvertype"  ).data(), "Klu" ) );
}

void SolverAmesos::setParameters()
{
    M_solver->SetParameters( M_trilinosParameterList );
}

// ===================================================
// Get Methods
// ===================================================
Int
SolverAmesos::NumIters()
{
    return 1;
}

Real
SolverAmesos::TrueResidual()
{
    return 0.;
}

// ===================================================
// Private Methods
// ===================================================
void SolverAmesos::createSolver( const std::string& solverType )
{
    Amesos factory;
    M_solver = factory.Create( solverType, M_problem );

    if ( M_solver == 0 )
    {
        if ( M_displayer.isLeader() )
        {
            std::cerr << std::endl  << std::endl;
            std::cerr << "SolverAmesos: Selected solver << " << solverType << " is not available. Bailing out." << std::endl;
        }

        // return ok not to break the test harness
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        exit( EXIT_SUCCESS );
    }
}

} // namespace LifeV

