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
    @brief LinearSolver

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 03-08-2011
 */

#ifndef _LINEARSOLVER_HPP
#define _LINEARSOLVER_HPP 1

#include <iomanip>


// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosEpetraOperator.h>
#include <BelosOutputManager.hpp>
#include <BelosMVOPTester.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosGCRODRSolMgr.hpp>
#include <BelosGmresPolySolMgr.hpp>
#include <BelosPCPGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosRCGSolMgr.hpp>
#include <BelosTFQMRSolMgr.hpp>

#include <Teuchos_ParameterList.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

#include <lifev/core/util/LifeDebug.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/operator/SolverOperator.hpp>
#include <lifev/core/operator/BelosOperator.hpp>
#include <lifev/core/operator/AztecooOperator.hpp>

namespace LifeV
{

//! LinearSolver - Class to wrap linear solver
/*!
  By default the solver is block gmres.

  @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
*/
class LinearSolver
{
public:

    //! @name Public Types
    //@{

    typedef Real                                                        value_Type;

    typedef LinearSolver                                                solver_Type;
    typedef Epetra_MultiVector                                          multiVector_Type;
    typedef boost::shared_ptr<multiVector_Type>                         multiVectorPtr_Type;
    typedef Epetra_Operator                                             operator_Type;
    typedef boost::shared_ptr<operator_Type>                            operatorPtr_Type;
    typedef Operators::SolverOperator                                   SolverOperator_Type;
    typedef boost::shared_ptr< SolverOperator_Type >                    SolverOperatorPtr_Type;
    typedef MatrixEpetra<Real>                                          matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                              matrixPtr_Type;
    typedef VectorEpetra                                                vector_Type;
    typedef boost::shared_ptr<VectorEpetra>                             vectorPtr_Type;
    typedef Preconditioner                                              preconditioner_Type;
    typedef boost::shared_ptr<preconditioner_Type>                      preconditionerPtr_Type;

    enum SolverType          { UndefinedSolver, Belos, AztecOO };

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    LinearSolver();

    //! Constructor
    /*!
      @param commPtr Communicator
     */
    LinearSolver( const boost::shared_ptr<Epetra_Comm> commPtr );

    //! Destructor
    ~LinearSolver();

    //@}

   //! @name Methods
    //@{

    //! Solves the system and returns the number of iterations.
    /*!
      The Matrix has already been passed by the method
      setMatrix or setOperator

      The preconditioner is build starting from the matrix baseMatrixForPreconditioner
      if it is set otherwise from the problem matrix.
      @param solutionPtr Vector to store the solution
      @return Number of iterations, M_maxIter+1 if solve failed.
     */
    Int solve( vectorPtr_Type solutionPtr );

    //! Compute the residual
    /*!
      @param solutionPtr Shared pointer on the solution of the system
      The method returns -1 if an error occurs
     */
    Real computeResidual( vectorPtr_Type solutionPtr );

    //! return the solver status
    std::string printStatus();

    //! Setup the preconditioner from a GetPot file
    /*!
      @param dataFile GetPot object which contains the data about the preconditioner
      @param section Section the GetPot structure where to find the informations about the preconditioner
     */
    void setPreconditionerFromGetPot( const GetPot& dataFile, const std::string& section );

    //! Builds the preconditioner starting from the matrix "baseMatrixForPreconditioner"
    /*!
      The preconditioner is build starting from the matrix baseMatrixForPreconditioner
      by the preconditioner object passed in by the method setPreconditioner
      @param  baseMatrixForPreconditioner Base matrix for the preconditioner construction
    */
    void buildPreconditioner();

    //! Reset the stored preconditioner
    /*!
      Note: This method only affects the LifeV::Preconditioner (i.e. not the Epetra_Operators
            used as preconditioner
     */
    void resetPreconditioner();

    //! Return true if preconditioner has been setted
    bool isPreconditionerSet() const;

    //! Reset the status for the state of convergence and loss of accuracy
    void resetStatus();

    //! Print informations about the solver
    void showMe( std::ostream& output = std::cout ) const;

    //! Setup the solver operator to be used
    void setupSolverOperator();

    //@}

    //! @name Set Method
    //@{

    //! Set the solver which should be used
    /*!
      @param solverOperatorType Type of solver manager
      The solver type can be chosen from one of the following:
      Aztecoo, Belos
     */
    void setSolverType( const SolverType& solverType );

    //! Method to set communicator for Displayer (for empty constructor)
    /*!
      @param commPtr Communicator for the displayer
     */
    void setCommunicator( const boost::shared_ptr<Epetra_Comm> commPtr );

    //! Method to set matrix from MatrixEpetra
    /*!
      @param matrixPtr Matrix of the system
     */
    void setOperator( matrixPtr_Type matrixPtr );

    //! Method to set a general linear operator (of class derived from Epetra_Operator) defining the linear system
    /*!
      @param operPtr Pointer to an operator for the system
     */
    void setOperator( operatorPtr_Type operPtr );

    //! Method to set the right hand side (rhs) of the linear system
    /*!
      @param rhsPtr right hand side of the system
     */
    void setRightHandSide( const vectorPtr_Type rhsPtr );

    //! Method to set an Preconditioner preconditioner
    /*!
      @param preconditionerPtr Preconditioner to be used to solve the system
     */
    void setPreconditioner( preconditionerPtr_Type preconditionerPtr );

    //! Method to set a general Epetra_Operator as preconditioner
    /*!
      @param preconditionerPtr  Preconditioner to be set of type Epetra_Operator
     */
    void setPreconditioner( operatorPtr_Type preconditionerPtr );

    //! Method to set a matrix on which the preconditioner should be created
    /*!
      @param baseMatrixPtr  matrix on which the preconditioner should be created
     */
    void setBaseMatrixForPreconditioner( matrixPtr_Type baseMatrixPtr );

    //! Method to setup the solver using Teuchos::ParameterList
    /*!
      @param list Teuchos::ParameterList object
      Note: The parameters are added to the existing one. Use resetParameters to clean the parameters list.
     */
    void setParameters( const Teuchos::ParameterList& list );

    //! Method to set a particular parameter
    /*!
      @param name Name of the parameter
      @param value Value of the parameter
      Note: The parameters are added to the existing one. Use resetParameters to clean the parameters list.
     */
    template<typename T>
    void setParameter( const std::string& name, T value );

    //! Method to reset the parameters list of the solver
    void resetParameters();

    //! Specify if the preconditioner should be reuse or not
    /*!
      @param reusePreconditioner If set to true, do not recompute the preconditioner
     */
    void setReusePreconditioner( const bool reusePreconditioner );

    //! Specify if the application should stop when problems occur in the iterations
    /*!
      @param enable If set to true, application will stop if problems occur
      Note: This option is useful for simulation on clusters. In particular if the
            system does not converge or if a loss of precision occurs time is saved
            by stoping the simulation
     */
    void setQuitOnFailure( const bool enable );

    //! Set the tolerance of the solver
    /*!
      @param tolerance Tolerance used by the solver
     */
    void setTolerance( const Real& tolerance );

    //@}

    //! @name Get Method
    //@{

    //! Return the total number of iterations
    Int numIterations() const;

    //! Return the recursive residual
    /*!
      The method returns -1 if an error occurs
     */
    Real recursiveResidual();

    //! Method to get a shared pointer to the preconditioner
    preconditionerPtr_Type& preconditioner();

    //! Return a Teuchos parameters list
    Teuchos::ParameterList& parametersList();

    //! Return a pointer on the Belos solver manager
    SolverOperatorPtr_Type solver();

    //! Return a shared pointer on the displayer
    boost::shared_ptr<Displayer> displayer();

    //! Returns the maximum of iterations tolerate to avoid recomputing the preconditioner
    Int maxItersForReuse() const;

    //! Returns if the preconditioner can be reused
    bool reusePreconditioner() const;

    //! Returns if the application should stop if a problem occurs
    bool quitOnFailure() const;

    //! Returns if the solver is in silent mode
    bool silent() const;

    //! Returns if the maximum number of iterations has been reached
    SolverOperator_Type::SolverOperatorStatusType hasReachedMaxNumIters() const;

    //! Returns if a loss of precision has been detected
    SolverOperator_Type::SolverOperatorStatusType isLossOfAccuracyDetected() const;

    //! Returns if the convergence has been achieved
    SolverOperator_Type::SolverOperatorStatusType hasConverged() const;

    //@}

private:

    //! @name Private Methods
    //@{

    //@}

    operatorPtr_Type             M_operator;
    matrixPtr_Type               M_matrix;
    matrixPtr_Type               M_baseMatrixForPreconditioner;
    vectorPtr_Type               M_rhs;

    preconditionerPtr_Type       M_preconditioner;
    operatorPtr_Type             M_preconditionerOperator;

    SolverType                   M_solverType;
    SolverOperatorPtr_Type       M_solverOperator;

    Teuchos::ParameterList       M_parameterList;
    boost::shared_ptr<Displayer> M_displayer;

    // LifeV features
    Int                          M_maxItersForReuse;
    bool                         M_reusePreconditioner;
    bool                         M_quitOnFailure;
    bool                         M_silent;

    // Status informations
    SolverOperator_Type::SolverOperatorStatusType M_lossOfPrecision;
    SolverOperator_Type::SolverOperatorStatusType M_maxNumItersReached;
    SolverOperator_Type::SolverOperatorStatusType M_converged;

    Real                        M_tolerance;

};

template<typename T>
void
LinearSolver::setParameter( const std::string& name, T value )
{
    M_parameterList.set( name, value );
}

} // namespace LifeV

#endif /* LINEARSOLVER_HPP */
