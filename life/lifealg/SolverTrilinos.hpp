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

#ifndef __SolverTrilinos_H
#define __SolverTrilinos_H 1

#include <iomanip>

#include <boost/shared_ptr.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <AztecOO_ConfigDefs.h>
#include <AztecOO.h>
#include <Teuchos_ParameterList.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifearray/EpetraVector.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraPreconditioner.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifecore/debug.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifecore/displayer.hpp>

namespace LifeV
{

//! SolverTrilinos - Class to wrap linear solver
/*!
  By default the solver is gmres and the preconditioner is ilu.

  @author Simone Deparis   <simone.deparis@epfl.ch>
  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
*/
class SolverTrilinos
{
public:

    //! @name Public Types
    //@{

    typedef Real                               value_type;

    typedef SolverTrilinos                     solver_type;

    typedef EpetraMatrix<Real>                 matrix_type;
    typedef EpetraVector                       vector_type;

    typedef EpetraPreconditioner               prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>   prec_type;
    typedef boost::shared_ptr<Epetra_Operator> comp_prec_type;
    typedef boost::shared_ptr<matrix_type>     matrix_ptrtype;
    typedef boost::shared_ptr<EpetraVector>    vector_ptrtype;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    SolverTrilinos();

    //! Constructor
    /*!
      @param comm Communicator
     */
    SolverTrilinos( const boost::shared_ptr<Epetra_Comm>& comm );

    //@}

   //! @name Methods
    //@{

    //! Solve the problem \f$ A x = b \f$.
    /*!
      A has been entered via \c setMatrix.
      @param solution Vector to store the solution
      @rhs rhs Right hand side of the problem
      @return Number of iterations, M_maxIter+1 if solve failed.
     */
    Int solve( vector_type& solution, const vector_type& rhs );

    //! Compute the residual
    /*!
      @param solution Solution of the system
      @param rhs Right hand side of the problem
     */
    Real computeResidual( vector_type& solution, vector_type& rhs );

    //! return the Aztec status
    std::string printStatus();

    //! Solves the system and returns the number of iterations.
    /*!
      The Matrix has already been passed by the method
      setMatrix or setOperator

      The preconditioner is build starting from the matrix baseMatrixForPreconditioner
      by the preconditioner object passed in by the method setPreconditioner
      @param  rhsFull Right hand side
      @param  solution vector to store the solution
      @param  baseMatrixForPreconditioner Base matrix for the preconditioner construction
      @return number of iterations. If negative, the solver did not converge,
               the preconditioner has been recomputed, and a second solution is tried
    */
    Int solveSystem( const vector_type& rhsFull,
                     vector_type&       solution,
                     matrix_ptrtype&    baseMatrixForPreconditioner );

    //! Solves the system and returns the number of iterations.
    /*!
      The Matrix has already been passed by the method
      setMatrix or setOperator
      @param  rhsFull Right hand side
      @param  solution Vector to store the solution
      @param  preconditionerPtr Pointer on a preconditioner to use (templated parameter, can derive from
                                EpetraPreconditioner class or from Epetra_Operator)
    */
    template <typename PrecPtrOperator>
    Int solveSystem(  const vector_type& rhsFull,
                      vector_type&       solution,
                      PrecPtrOperator    preconditionerPtr );

    //! Setup the preconditioner
    /*!
      @param dataFile GetPot object which contains the data about the preconditioner
      @param section Section the GetPot structure where to find the informations about the preconditioner
     */
    void setUpPrec( const GetPot& dataFile, const std::string& section );

    //! Builds the preconditioner starting from the matrix "baseMatrixForPreconditioner"
    /*!
      The preconditioner is build starting from the matrix baseMatrixForPreconditioner
      by the preconditioner object passed in by the method setPreconditioner
      @param  baseMatrixForPreconditioner Base matrix for the preconditioner construction
    */
    void buildPreconditioner( matrix_ptrtype& baseMatrixForPreconditioner );

    //! Delete the stored preconditioner
    void precReset();

    //! Return true if preconditioner has been setted
    bool isPrecSet() const;

    //! Print informations about the solver
    void showMe( std::ostream& output = std::cout ) const;

    //@}

    //! @name Set Method
    //@{

    //! Method to set communicator for Displayer (for empty constructor)
    /*!
      @param comm Communicator for the displayer
     */
    void setCommunicator( const boost::shared_ptr<Epetra_Comm>& comm );

    //! Method to set matrix from EpetraMatrix
    /*!
      @param matrix Matrix of the system
     */
    void setMatrix( matrix_type& matrix );

    //! Method to set a general linear operator (of class derived from Epetra_Operator) defining the linear system
    /*!
      @param oper Operator for the system
     */
    void setOperator( Epetra_Operator& oper );

    //! Method to set an EpetraPreconditioner preconditioner
    /*!
      @param preconditioner EpetraPreconditioner to be used to solve the system
     */
    void setPreconditioner( prec_type& preconditioner );

    //! Method to set a general Epetra_Operator as preconditioner
    /*!
      @param preconditioner  Preconditioner to be set of type Epetra_Operator
     */
    void setPreconditioner( comp_prec_type& preconditioner );

    //! Method to setup the solver using GetPot
    /*!
      @param dataFile GetPot object which contains the data about the solver
     */
    void setDataFromGetPot( const GetPot& dataFile, const std::string& section );

    //! Set the current parameters with the internal parameters list
    /*!
      Note: The parameter list is set using "setDataFromGetPot".
      @param cerrWarningIfUnused If true the solver return warning if some parameters are unused
     */
    void setParameters( bool cerrWarningIfUnused = false );

    //! Set the tolerance and the maximum number of iteration
    /*!
      @param tolerance Tolerance for the solver
      @param maxIter Maximum number of iteration
     */
    void setTolMaxIteration( const Real tolerance, const Int maxIter = -1 );

    //! Specify if the preconditioner should be reuse or not
    /*!
      @param reusePreconditioner If set to true, do not recompute the preconditioner
     */
    void setReusePreconditioner( const bool reusePreconditioner );

    //! Return the displayer
    boost::shared_ptr<Displayer> displayer();

    //@}

    //! @name Get Method
    //@{

    //! Return the total number of iterations
    Int NumIters() const;

    //! Return the maximum total number of iterations
    Int MaxIter() const;

    //! Return the true residual
    Real TrueResidual();

    //! Method to get a shared pointer to the preconditioner (of type derived from EpetraPreconditioner)*/
    prec_type& getPrec();

    //! Return the AztecStatus
    void getAztecStatus( Real status[AZ_STATUS_SIZE] );

    //! Return a Teuchos parameters list
    Teuchos::ParameterList& getParameterList();

    //! Return a reference on the AztecOO solver
    AztecOO& getSolver();

    //@}

private:

    matrix_type::matrix_ptrtype  M_matrix;
    prec_type                    M_preconditioner;

    AztecOO                      M_solver;

    Teuchos::ParameterList       M_TrilinosParameterList;
    boost::shared_ptr<Displayer> M_displayer;

    Real                         M_tolerance;
    Int                          M_maxIter;
    Int                          M_maxIterForReuse;
    bool                         M_reusePreconditioner;
};

template <typename PrecPtrOperator>
Int SolverTrilinos::solveSystem( const vector_type&  rhsFull,
                                 vector_type&        solution,
                                 PrecPtrOperator     preconditioner )

{
    M_displayer->leaderPrint("      AztecOO solving system ...         ");
    setPreconditioner(preconditioner);

    LifeChrono chrono;
    chrono.start();
    Int numIter = solve( solution, rhsFull );
    chrono.stop();
    M_displayer->leaderPrintMax( "      done in " , chrono.diff() );

    // If we use the "none" as output setting, we display just a summary
    if ( M_TrilinosParameterList.get( "output", "all" ) == "none" )
    {
        M_displayer->leaderPrint( "      Iterations number:                       ", M_solver.NumIters(), "\n" );
        M_displayer->leaderPrint( "      Scaled residual:                         ", M_solver.ScaledResidual(), "\n" );
    }

    if ( numIter >= M_maxIter )
        numIter = -numIter;

    return numIter;
}


} // namespace LifeV

#endif /* __SolverTrilinos_H */
