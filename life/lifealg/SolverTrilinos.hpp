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

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <AztecOO_config.h>
#include <AztecOO.h>
#include <Teuchos_ParameterList.hpp>

#include <life/lifearray/EpetraVector.hpp>
#include <life/lifearray/EpetraMatrix.hpp>

#include <life/lifealg/EpetraPreconditioner.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>

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
    typedef boost::shared_ptr<Epetra_Operator> comp_prec_type;
    typedef boost::shared_ptr<matrix_type>   matrix_ptrtype;
    typedef boost::shared_ptr<EpetraVector>  vector_ptrtype;


    //! @name Constructors & Destructor
    //@{

    SolverTrilinos();
    SolverTrilinos( const boost::shared_ptr<Epetra_Comm> );

    //@}


    //! @name Get Method
    //@{

    //! Total number of iterations
    int NumIters() const;

    //! Maximum Total number of iterations
    int MaxIter() const;

    //! True Residual
    double TrueResidual();

    /** Method to get a shared pointer to the preconditioner (of type derived from EpetraPreconditioner)*/
    prec_type& getPrec();

    void getAztecStatus( double status[AZ_STATUS_SIZE]);

    Teuchos::ParameterList& getParameterList();

    AztecOO& getSolver();

    //@}


    //! @name Set Method
    //@{

    //! Method to set communicator for Displayer (for empty constructor)
    void SetCommunicator( const boost::shared_ptr<Epetra_Comm> comm);

    //! Method to set matrix from EpetraMatrix
    void setMatrix(matrix_type& m);

    /** Method to set a general linear operator (of class derived from Epetra_Operator) defining the linear system*/
    void setOperator(Epetra_Operator& op);

    //! Method to set an EpetraPreconditioner preconditioner
    void setPreconditioner( prec_type& _prec );

    /** Method to set a general Epetra_Operator as preconditioner*/
    void setPreconditioner( comp_prec_type& _prec );

    void setDataFromGetPot( const GetPot& dfile, const std::string& section );

    void SetParameters( bool cerr_warning_if_unused = false );

    void setTolMaxiter( const double tol, const int maxiter = -1 );

    //! if set to true,  do not recompute the preconditioner
    void setReusePreconditioner( const bool reuse );

    boost::shared_ptr<Displayer> displayer(){return M_displayer;}
    //@}


    //! @name Methods
    //@{

    //! solve
    /*! Solve the problem \f$ A x = b \f$. \c A has been entered via \c setMatrix.
     *  return the number of iterations, M_maxIter+1 if solve failed.
     * \param algorithm - MS_Algorithm
     */
    int solve( vector_type& x, const vector_type& b );

    double computeResidual( vector_type& __X, vector_type& __B );

    // return the Aztec status
    std::string printStatus();

    //! Solves the system and returns the number of iterations.
    /*! The Matrix has already been passed by the method
        setMatrix or setOperator
        The preconditioner is build starting from the matrix baseMatrixForPreconditioner
        by the preconditioner object passed in by the method setPreconditioner
        @param  rhsFull   right hand side
        @param  sol       solution
        @param  baseMatrixForPreconditioner base matrix for the preconditioner construction

        returns number of iterations. If negative, the solver did not converge,
        the preconditioner has been recomputed, and a second solution is tried
    */
    int solveSystem(  const vector_type& rhsFull,
                      vector_type&       sol,
                      matrix_ptrtype&    baseMatrixForPreconditioner );

    //! Solves the system and returns the number of iterations.
    /*! The Matrix has already been passed by the method
        setMatrix or setOperator
        @param  rhsFull   right hand side
        @param  sol       solution
        @param  prec      preconditioner to use (templated parameter, can derive from
EpetraPreconditioner class or from Epetra_Operator)
    */
    template <typename PrecPtrOperator>
    int solveSystem(  const vector_type& rhsFull,
                      vector_type&       sol,
                      PrecPtrOperator         prec );

    void setUpPrec(const GetPot& dataFile,  const std::string& section);

    //! builds the preconditioner starting from the matrix baseMatrixForPreconditioner
    /*! The preconditioner is build starting from the matrix baseMatrixForPreconditioner
        by the preconditioner object passed in by the method setPreconditioner
        @param  baseMatrixForPreconditioner base matrix for the preconditioner construction
    */
    void buildPreconditioner( matrix_ptrtype& baseMatrixForPreconditioner);


    //@}


    // Return if preconditioner has been setted
    bool isPrecSet() const;

    void precReset() { M_prec->precReset(); }

private:

    matrix_type::matrix_ptrtype M_matrix;
    prec_type                   M_prec;

    AztecOO                     M_solver;

    Teuchos::ParameterList      M_TrilinosParameterList;
    boost::shared_ptr<Displayer> M_displayer;

    double                      M_tol;
    int                         M_maxIter;
    int                         M_maxIterForReuse;
    bool                        M_reusePreconditioner;
};

template <typename PrecPtrOperator>
int SolverTrilinos::solveSystem( const vector_type&  rhsFull,
                                 vector_type&        sol,
                                 PrecPtrOperator          prec )

{
    M_displayer->leaderPrint("      Solving system ...                   \n");

    setPreconditioner(prec);

    Chrono chrono;
    chrono.start();
    int numIter = solve( sol, rhsFull );
    chrono.stop();
    M_displayer->leaderPrintMax( "       ... done in " , chrono.diff() );

    if (numIter >= M_maxIter)
        numIter = -numIter;

    return numIter;
}


} // namespace LifeV

#endif /* __SolverTrilinos_H */
