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
    @brief SolverBelos

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 03-08-2011
 */

#ifndef SOLVERBELOS_H
#define SOLVERBELOS_H 1

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

#include <life/lifearray/VectorEpetra.hpp>
#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifealg/Preconditioner.hpp>
#include <life/lifecore/LifeDebug.hpp>
#include <life/lifefilters/GetPot.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifecore/Displayer.hpp>

namespace LifeV
{

//! SolverBelos - Class to wrap linear solver
/*!
  By default the solver is block gmres.

  @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
*/
class SolverBelos
{
public:

    //! @name Public Types
    //@{

    typedef Real                                                        value_Type;

    typedef SolverBelos                                                 solver_Type;
    typedef Epetra_MultiVector                                          multiVector_Type;
    typedef Epetra_Operator                                             operator_Type;
    typedef boost::shared_ptr<operator_Type>                            operatorPtr_Type;
    typedef Belos::SolverManager<Real,multiVector_Type,operator_Type>   SolverManager_Type;
    typedef RCP< SolverManager_Type >                                   SolverManagerPtr_Type;
    typedef Belos::LinearProblem<double,multiVector_Type,operator_Type> LinearProblem_Type;
    typedef RCP< LinearProblem_Type >                                   LinearProblemPtr_Type;

    typedef MatrixEpetra<Real>                                          matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                              matrixPtr_Type;
    typedef VectorEpetra                                                vector_Type;
    typedef boost::shared_ptr<VectorEpetra>                             vectorPtr_Type;
    typedef Preconditioner                                              preconditioner_Type;
    typedef boost::shared_ptr<preconditioner_Type>                      preconditionerPtr_Type;



    enum PrecApplicationType { LeftPreconditioner, RightPreconditioner };
    enum SolverManagerType { BlockCG, PseudoBlockCG, RCG,
                             BlockGmres, PseudoBlockGmres, BlockFGmres, PseudoBlockFGmres, GmresPoly,
                             GCRODR, PCPG, TFQMR };

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    SolverBelos();

    //! Constructor
    /*!
      @param comm Communicator
     */
    SolverBelos( const boost::shared_ptr<Epetra_Comm>& comm );

    //! Destructor
    ~SolverBelos();

    //@}

   //! @name Methods
    //@{

    //! Solves the system and returns the number of iterations.
    /*!
      The Matrix has already been passed by the method
      setMatrix or setOperator

      The preconditioner is build starting from the matrix baseMatrixForPreconditioner
      if it is set otherwise from the problem matrix.
      @param solution Vector to store the solution
      @return Number of iterations, M_maxIter+1 if solve failed.
     */
    Int solve( vector_Type& solution );

    //! Solves the system and returns the number of iterations.
    /*!
      The Matrix has already been passed by the method
      setMatrix or setOperator

      The preconditioner is build starting from the matrix baseMatrixForPreconditioner
      if it is set otherwise from the problem matrix.
      @param solution Vector to store the solution
      @return Number of iterations, M_maxIter+1 if solve failed.
     */
    Int solve( multiVector_Type& solution );

    //! Compute the residual
    /*!
      @param solution Solution of the system
      The method returns -1 if an error occurs
     */
    Real computeResidual( vector_Type& solution );

    //! return the solver status
    std::string printStatus();

    //! Setup the preconditioner from a GetPot file
    /*!
      @param dataFile GetPot object which contains the data about the preconditioner
      @param section Section the GetPot structure where to find the informations about the preconditioner
     */
    void setPreconditionerFromGetPot( const GetPot& dataFile, const std::string& section, PrecApplicationType precType = RightPreconditioner );

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

    //! Print informations about the solver
    void showMe( std::ostream& output = std::cout ) const;

    //@}

    //! @name Set Method
    //@{

    //! Set the solver manager which should be used by Belos
    /*!
      @param solverManager type of solver manager
      The solver manager can be chosen from one of the following:
      BlockCG, PseudoBlockCG, RCG, BlockGmres, PseudoBlockGmres,
      BlockFGmres, PseudoBlockFGmres, GCRODR, GmresPoly,PCPG, TFQMR.
     */
    void setSolverManager( const SolverManagerType& solverManager );

    //! Method to set communicator for Displayer (for empty constructor)
    /*!
      @param comm Communicator for the displayer
     */
    void setCommunicator( const boost::shared_ptr<Epetra_Comm>& comm );

    //! Method to set matrix from MatrixEpetra
    /*!
      @param matrix Matrix of the system
     */
    void setMatrix( matrixPtr_Type& matrix );

    //! Method to set a general linear operator (of class derived from Epetra_Operator) defining the linear system
    /*!
      @param oper Operator for the system
     */
    void setOperator( Epetra_Operator& oper );

    //! Method to set the right hand side (rhs) of the linear system
    /*!
      @param rhs right hand side of the system
     */
    void setRightHandSide( const vector_Type& rhs );

    //! Method to set the right hand side (rhs) of the linear system
    /*!
      @param rhs right hand side of the system
     */
    void setRightHandSide( const multiVector_Type& rhs );

    //! Method to set an Preconditioner preconditioner
    /*!
      @param preconditioner Preconditioner to be used to solve the system
     */
    void setPreconditioner( preconditionerPtr_Type& preconditioner, PrecApplicationType precType = RightPreconditioner );

    //! Method to set a general Epetra_Operator as preconditioner
    /*!
      @param preconditioner  Preconditioner to be set of type Epetra_Operator
     */
    void setPreconditioner( operatorPtr_Type& preconditioner, PrecApplicationType precType = RightPreconditioner );

    //! Method to setup the solver using GetPot
    /*!
      @param dataFile GetPot object which contains the data about the solver
      Note: The parameters are added to the existing one. Use resetParameters to clean the parameters list.
     */
    void setParameters( const GetPot& dataFile, const std::string& section );

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

    //@}

    //! @name Get Method
    //@{

    //! Return the total number of iterations
    Int numIterations() const;

    //! Return the true residual
    /*!
      The method returns -1 if an error occurs
     */
    Real trueResidual();

    //! Method to get a shared pointer to the preconditioner (of type derived from Preconditioner)*/
    preconditionerPtr_Type& preconditioner( PrecApplicationType precType = RightPreconditioner );

    //! Return a Teuchos parameters list
    Teuchos::ParameterList& parametersList();

    //! Return a pointer on the Belos solver manager
    SolverManagerPtr_Type solver();

    //! Return a shared pointer on the displayer
    boost::shared_ptr<Displayer> displayer();

    //@}

private:

    //! @name Private Methods
    //@{

    //! Setup the solver manager to be used
    void setupSolverManager();

    //@}

    matrixPtr_Type               M_matrix;
    matrixPtr_Type               M_baseMatrixForPreconditioner;
    preconditionerPtr_Type       M_leftPreconditioner;
    preconditionerPtr_Type       M_rightPreconditioner;

    SolverManagerType            M_solverManagerType;
    SolverManagerPtr_Type        M_solverManager;
    LinearProblemPtr_Type        M_problem;

    Teuchos::ParameterList       M_parameterList;
    boost::shared_ptr<Displayer> M_displayer;

    // LifeV features for Belos
    Int                          M_maxItersForReuse;
    bool                         M_reusePreconditioner;
    bool                         M_quitOnFailure;
    bool                         M_silent;

    // Status information
    bool                         M_lossOfPrecision;
    bool                         M_maxNumItersReached;

};

template<typename T>
void
SolverBelos::setParameter( const std::string& name, T value )
{
    M_parameterList.set( name, value );
}

} // namespace LifeV

#endif /* SOLVERBELOS_H */
