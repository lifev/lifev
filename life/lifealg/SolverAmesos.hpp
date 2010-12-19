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

#ifndef _SolverAmesos_H
#define _SolverAmesos_H

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Amesos.h>
#include <Amesos_BaseSolver.h>
#include <Teuchos_ParameterList.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifearray/EpetraVector.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifecore/displayer.hpp>

class GetPot;

namespace LifeV
{

//! SolverAmesos - Class to wrap linear solver
/*!
  @author Simone Deparis   <simone.deparis@epfl.ch>
  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
*/
class SolverAmesos
{
public:

    //! @name Public Types
    //@{

    typedef Real                             value_type;

    typedef Displayer::commPtr_Type          commPtr_Type;

    typedef SolverAmesos                     solver_type;

    typedef EpetraMatrix<Real>               matrix_type;
    typedef EpetraVector                     vector_type;

    typedef void                             prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type> prec_type;
    typedef boost::shared_ptr<matrix_type>   matrix_ptrtype;
    typedef boost::shared_ptr<EpetraVector>  vector_ptrtype;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Default constructor
    /*!
     * @param comm The communicator.
     */
    SolverAmesos( const commPtr_Type& comm );

    //! Destructor
    ~SolverAmesos() {}

    //@}


    //! @name Methods
    //@{

    //! Compute the residual
    /*!
      @param solution Solution vector
      @param rhs Right hand side of the system
     */
    Real computeResidual( const vector_type& solution, const vector_type& rhs );

    //! Solves the system and returns the number of iterations.
    /*!
      returns number of iterations. If negative, the solver did not converge,
      the preconditionar has been recomputed, and a second solution is tried
      @param  rhsFull Right hand side vector
      @param  solution Vector to store the solution
    */
    Int solveSystem( vector_type&    rhsFull,
                     vector_type&    solution,
                     matrix_ptrtype& /* unused */ );

    //! Display status of the solver
    void printStatus();

    //! Return true if the preconditioner is set
    /*!
      Note: This method always return true!
    */
    bool isPrecSet() const;

    //! Delete the stored preconditioner
    /*!
      Note: This method is empty
     */
    void precReset();

    //! Setup the preconditioner
    /*!
      Note: This method is empty
      @param dataFile GetPot object which contains the data about the preconditioner
      @param section Section the GetPot structure where to find the informations about the preconditioner
     */
    void setUpPrec( const GetPot& dataFile,  const std::string& section );

    //! Specify if the preconditioner should be reuse or not
    /*!
      Note: This method is empty
      @param reusePreconditioner If set to true, do not recompute the preconditioner
     */
    void setReusePreconditioner( const bool& /*reusePreconditioner*/ );

    //! Print informations about the solver
    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name  Set Methods
    //@{

    //! Set matrix from EpetraMatrix
    /*!
      @param matrix Matrix of the system
     */
    void setMatrix( const matrix_type& matrix );

    //! Method to set a general linear operator (of class derived from Epetra_Operator) defining the linear system
    /*!
      @param oper Operator for the system
     */
    void setOperator( const Epetra_Operator& oper );

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
    void setParameters();

    //! Set the tolerance and the maximum number of iteration
    /*!
      @param tolerance Tolerance for the solver
      @param maxIter Maximum number of iteration
     */
    void setTolMaxIteration( const Real tolerance, const Int maxIter = -1 );

    //@}


    //! @name Get Methods
    //@{

    //! Return the total number of iterations
    Int NumIters();

    //! Return the true residual
    Real TrueResidual();

    //@}

private:

    //! @name Private Methods
    //@{

    //! Create a solver using a factory
    /*!
      @param solverType String containing the name of the solver
     */
    void createSolver( const std::string& solverType );

    //@}

    matrix_type::matrix_ptrtype M_matrix;

    Epetra_LinearProblem        M_problem;

    Amesos_BaseSolver*          M_solver;

    Teuchos::ParameterList      M_trilinosParameterList;

    Displayer                   M_displayer;
};

} // namespace LifeV

#endif /* _SolverAmesos_H */
