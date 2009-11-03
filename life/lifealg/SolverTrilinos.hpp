/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Simone Deparis   <simone.deparis@epfl.ch>
            Gilles Fourestey <gilles.fourestey@epfl.ch>
      Date: 2006-11-08

 Copyright (C) 2006 EPFL

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
   \author Simone Deparis   <simone.deparis@epfl.ch>
   \author Gilles Fourestey <gilles.fourestey@epfl.ch>
   \date 2006-11-08
*/

#ifndef __SolverTrilinos_H
#define __SolverTrilinos_H 1

#include "Epetra_ConfigDefs.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "AztecOO_config.h"
#include "AztecOO.h"
#include "Teuchos_ParameterList.hpp"

#include "life/lifearray/EpetraVector.hpp"
#include "life/lifearray/EpetraMatrix.hpp"

#include "life/lifealg/EpetraPreconditioner.hpp"
#include "life/lifealg/IfpackPreconditioner.hpp"

#include <life/lifecore/debug.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifecore/displayer.hpp>

#include <boost/shared_ptr.hpp>

#include <iomanip>

namespace LifeV
{
// namespace Epetra
// {
/*!
  \class SolverTrilinos
  \brief wrap trilinos linear solvers

  By default the solver is gmres and the preconditioner is ilu.

  \author Simone Deparis   <simone.deparis@epfl.ch>
  \author Gilles Fourestey <gilles.fourestey@epfl.ch>
  @see
*/

class SolverTrilinos
{
public:

    typedef double                           value_type;

    typedef SolverTrilinos                   solver_type;

    typedef EpetraMatrix<double>             matrix_type;
    typedef EpetraVector                     vector_type;

    typedef EpetraPreconditioner             prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type> prec_type;
    typedef boost::shared_ptr<matrix_type>   matrix_ptrtype;
    typedef boost::shared_ptr<EpetraVector>  vector_ptrtype;


    //! @name Constructors & Destructor
    //@{

    SolverTrilinos();
    SolverTrilinos( const Epetra_Comm& );

    //@}


    //! @name Get Method
    //@{

    //! Total number of iterations
    int    NumIters();

    //! True Residual
    double TrueResidual();

    // Return if preconditioner has been setted
    bool precSet() const
    {
        return ( M_prec.get() !=0 && M_prec->getPrec() != 0 );
    }

    prec_type& getPrec()
    {
        return M_prec;
    }

    void getAztecStatus( double status[AZ_STATUS_SIZE])
    {
        M_solver.GetAllAztecStatus( status );
    }

    //@}


    //! @name Set Method
    //@{

    //! Set communicator for Displayer (for empty constructor)
    void SetCommunicator( const Epetra_Comm& comm);

    //! set matrix from EpetraMatrix
    void setMatrix(matrix_type& m);

    void setOperator(Epetra_Operator& op);

    //! set Epetra_Operator preconditioner
    void setPreconditioner( prec_type _prec );

    void setPrec( prec_raw_type* prec );

    void setDataFromGetPot( const GetPot& dfile, const std::string& section );
    void SetParameters( bool cerr_warning_if_unused = false );

    void useGMRES(const int restart=300);
    void useCG();
    void useBICGSTAB();
    void useJacobyPrec();
    void useNoPrec();
    void useDDPrec();
    void useNeumannPrec();
    void setReuse();
    void setTolMaxiter( const double tol, const int maxiter = -1 );

    enum VerboseLevel { NONE, SUMMARY, LAST }; //! describes the verobisty level
    void SetVerbose( const VerboseLevel verb=NONE );

    //@}


    //! @name Methods
    //@{

    //! solve
    /*! Solve the problem \f$ A x = b \f$. \c A has been entered via \c setMatrix.
     *  return the number of iterations, M_maxIter+1 if solve failed.
     * \param algorithm - MS_Algorithm
     */
    int solve( vector_type& x, vector_type& b );

    double computeResidual( vector_type& __X, vector_type& __B );

    // return the Aztec status
    std::string printStatus();

    //@}

    /** Solves the system and returns the number of iterations.
        @param  matrFull,
        @param  rhsFull,
        @param  sol,
        @param  prec,
        @param  reuse,
        @param  retry = true

        returns number of iterations. If negative, the solver did not converge,
        the preconditionar has been recomputed, and a second solution is tried
    */
    template<typename PrecType>
    int solveSystem(  vector_type&     rhsFull,
                      vector_type&     sol,
                      PrecType&        prec,
                      bool const       reuse,
                      bool const       retry=true  );

//     int solveSystem(  vector_type&     rhsFull,
//                       vector_type&     sol,
//                       matrix_ptrtype&  prec,
//                       bool const       reuse,
//                       bool const       retry=true);

    void setUpPrec(const GetPot& dataFile,  const std::string& section);

    void buildPreconditioner( matrix_ptrtype& prec);

    template<typename PrecType>
    void buildPreconditioner( PrecType& prec) { M_prec = prec; }

    void precReset() { M_prec->precReset(); }

    void precReset(prec_type& prec);

    template<typename PrecType>
    void precReset(PrecType& prec) {}

private:

    matrix_type::matrix_ptrtype M_matrix;
    prec_type                   M_prec;

    AztecOO                     M_solver;

    Teuchos::ParameterList      M_TrilinosParameterList;
    Displayer                   M_Displayer;

    int                         M_maxIter;
    double                      M_tol;
    int                         M_maxIterSolver;
    int                         M_maxIterForReuse;
};

template<typename PrecType>
int SolverTrilinos::solveSystem(  vector_type&      rhsFull,
                                  vector_type&      sol,
                                  PrecType&        prec,
                                  bool const        reuse,
                                  bool const        retry)

{
    M_Displayer.leaderPrint("      Setting up the solver ...                \n");

    if ( !M_prec->set() || !reuse  )
    {
        buildPreconditioner(prec);
        setPreconditioner(M_prec);
    }
    else
    {
        M_Displayer.leaderPrint("      Reusing  precond ...                 \n");
    }

    M_Displayer.leaderPrint("      Solving system ...                   \n");

    Chrono chrono;
    chrono.start();
    int numIter = solve( sol, rhsFull );
    chrono.stop();
    M_Displayer.leaderPrintMax( "       ... done in " , chrono.diff() );

    // If we do not want to retry, return now.
    // otherwise rebuild the preconditioner and solve again:
    if ( numIter >= M_maxIterSolver && retry )
    {
        chrono.start();

        M_Displayer.leaderPrint("     Iterative solver failed, numiter = " , numIter);
        M_Displayer.leaderPrint("     maxIterSolver = " , M_maxIterSolver );
        M_Displayer.leaderPrint("     retrying:          ");

        if (prec.get())
        {
            precReset(M_prec);
            buildPreconditioner(prec);
            setPreconditioner(M_prec);
        }
        chrono.stop();
        M_Displayer.leaderPrintMax( "done in " , chrono.diff() );
        // Solving again, but only once (retry = false)
        numIter = solveSystem( rhsFull, sol, prec, reuse, false );

        if (numIter >= M_maxIterSolver)
            M_Displayer.leaderPrint(" ERROR: Iterative solver failed again.\n");
        return -numIter;
    }
    return numIter;

//     Chrono chrono;
//     setPreconditioner(prec);
//     chrono.start();
//     int numIter = solve(sol, rhsFull);
//     chrono.stop();
//     M_Displayer.leaderPrintMax( "       ... done in " , chrono.diff() );
//     return numiter;
}

} // namespace LifeV

#endif /* __SolverTrilinos_H */

/*
//! set the level of recursion in the solution of a linear system
void setRecursionLevel( int newLevel )
{
    M_options[ AZ_recursion_level ] = newLevel;
}
*/
