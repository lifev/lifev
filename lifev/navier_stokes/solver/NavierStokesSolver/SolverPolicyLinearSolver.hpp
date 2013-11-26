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
    @file SolverPolicyLinearSolver class
    @brief This class contains all the informations necessary
           to solve a linear system

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 07-12-2012
 */

#ifndef SOLVERPOLICYLINEARSOLVER_HPP
#define SOLVERPOLICYLINEARSOLVER_HPP

#include <iostream>
#include <boost/shared_ptr.hpp>


#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>


#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>


namespace LifeV
{

struct SolverPolicyLinearSolver
{
public:
    typedef MatrixEpetra<Real>                       matrix_Type;
    typedef boost::shared_ptr<matrix_Type>           matrixPtr_Type;
    typedef VectorEpetra                             vector_Type;
    typedef boost::shared_ptr<VectorEpetra>          vectorPtr_Type;
    typedef Epetra_Comm                              comm_Type;
    typedef boost::shared_ptr<comm_Type>             commPtr_Type;
    typedef LinearSolver                             solver_Type;
    typedef boost::shared_ptr< solver_Type >         solverPtr_Type;
    typedef Preconditioner                           preconditioner_Type;
    typedef boost::shared_ptr<preconditioner_Type>   preconditionerPtr_Type;

    //! Method to set a preconditioner
    /*!
        @param preconditionerPtr Preconditioner to be used to solve the system
     */
    void setPreconditioner ( preconditionerPtr_Type preconditionerPtr );

    //! Method to get a preconditioner
    preconditionerPtr_Type preconditioner();

protected:
    void initSolver ( Teuchos::ParameterList& list );

    int solve ( matrixPtr_Type systemMatrix,
                vectorPtr_Type rhs,
                vectorPtr_Type solution );
    solverPtr_Type M_solver;

    virtual Displayer displayer() = 0;
    virtual commPtr_Type comm() = 0;
};

} // namespace LifeV

#endif /* EXPORTERPOLICYHDF5_HPP */
