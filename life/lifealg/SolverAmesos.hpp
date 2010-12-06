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

#include <Amesos.h>
#include <Amesos_BaseSolver.h>

#include <Teuchos_ParameterList.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <life/lifearray/EpetraMatrix.hpp>

#include <life/lifecore/displayer.hpp>

class GetPot;

namespace LifeV
{

class SolverAmesos
{
public:

    //! @name Public Types
    //@{

    typedef Real                             value_type;

    typedef Displayer::comm_PtrType          comm_PtrType;

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
     * @param comm the communicator.
     */
    SolverAmesos( const comm_PtrType& comm );

    //! Destructor
    ~SolverAmesos() {}

    //@}


    //! @name Methods
    //@{

    Real computeResidual( const vector_type& x, const vector_type& rhs );

    /** Solves the system and returns the number of iterations.
        @param  matrFull,
        @param  rhsFull,
        @param  sol,
        @param  prec
        returns number of iterations. If negative, the solver did not converge,
        the preconditionar has been recomputed, and a second solution is tried
    */
    Int solveSystem( vector_type&    rhsFull,
                     vector_type&    sol,
                     matrix_ptrtype& /* unused */);

    //Display status
    void printStatus();

    bool isPrecSet() const;

    void precReset();

    void setUpPrec( const GetPot& dataFile,  const std::string& section );

    void setReusePreconditioner( const bool& /*reuse*/ );

    //@}


    //! @name  Set Methods
    //@{

    //! Set matrix from EpetraMatrix
    void setMatrix( const matrix_type& matrix );

    void setOperator( const Epetra_Operator& op );

    void setDataFromGetPot( const GetPot& dfile, const std::string& section );

    void setParameters();

    void setTolMaxiter( const Real tol, const Int maxiter = -1 );

    //@}


    //! @name Get Methods
    //@{

    //! Total number of iterations
    Int NumIters();

    //! True Residual
    Real TrueResidual();

    //@}

private:

    //! @name Private Methods
    //@{

    void createSolver( const std::string& solverType );

    //@}

    matrix_type::matrix_ptrtype M_matrix;

    Epetra_LinearProblem        M_problem;

    Amesos_BaseSolver*          M_solver;

    Teuchos::ParameterList      M_TrilinosParameterList;

    Displayer                   M_displayer;
};

} // namespace LifeV

#endif /* _SolverAmesos_H */
