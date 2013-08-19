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
    @file TimeIterationNonlinear class
    @brief This class contains all the informations necessary
           to perform a time iteration

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 06-12-2012
 */

#ifndef TIMEITERATIONPOLICYNONLINEAR_HPP
#define TIMEITERATIONPOLICYNONLINEAR_HPP

#include <iostream>
#include <string>
#include <boost/shared_ptr.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/SolverPolicyLinearSolver.hpp>


namespace LifeV
{

template< class AssemblyPolicy, class SolverPolicy = SolverPolicyLinearSolver >
struct TimeIterationPolicyNonlinear : private AssemblyPolicy, public virtual SolverPolicy
{
protected:
    typedef MatrixEpetra<Real>                       matrix_Type;
    typedef boost::shared_ptr<matrix_Type>           matrixPtr_Type;
    typedef VectorEpetra                             vector_Type;
    typedef boost::shared_ptr<VectorEpetra>          vectorPtr_Type;
    typedef RegionMesh<LinearTetra>                  mesh_Type;
    typedef MeshPartitioner< mesh_Type >             meshPartitioner_Type;
    typedef MapEpetra                                map_Type;
    typedef boost::shared_ptr<map_Type>              mapPtr_Type;
    typedef FESpace< mesh_Type, map_Type >           fespace_Type;
    typedef boost::shared_ptr< fespace_Type >        fespacePtr_Type;
    typedef TimeAdvanceBDF<vector_Type>              bdf_Type;
    typedef boost::shared_ptr< bdf_Type >            bdfPtr_Type;
    typedef BCHandler                                bcContainer_Type;
    typedef boost::shared_ptr<bcContainer_Type>      bcContainerPtr_Type;

    enum { BDFOrder = AssemblyPolicy::BDFOrder };

    void initTimeIteration ( Teuchos::ParameterList& list );

    void iterate ( vectorPtr_Type solution,
                   bcContainerPtr_Type bchandler,
                   const Real& currentTime );

    // Parameters
    bool M_computeResidual;
    Real M_nonLinearTolerance;

    matrixPtr_Type          M_systemMatrix;
    mapPtr_Type             M_solutionMap;
    vectorPtr_Type          M_rhs;

    virtual Displayer displayer() = 0;
    virtual fespacePtr_Type uFESpace() const = 0;
    virtual fespacePtr_Type pFESpace() const = 0;
};

template< class AssemblyPolicy, class SolverPolicy >
void
TimeIterationPolicyNonlinear<AssemblyPolicy, SolverPolicy>::
initTimeIteration ( Teuchos::ParameterList& list )
{
    // Parameters
    M_computeResidual = list.get ( "Compute exact residual", false );
    M_nonLinearTolerance = list.get ( "Nonlinear tolerance", 1e-6 );

    // Initialization
    M_solutionMap.reset ( new map_Type ( uFESpace()->map() + pFESpace()->map() ) );

    // Init the assembler
    Teuchos::ParameterList assemblyList = list.sublist ( "Assembly: Parameter list" );
    AssemblyPolicy::initAssembly ( assemblyList );

    // Init the solver
    Teuchos::ParameterList solverList = list.sublist ( "Solver: Parameter list" );
    SolverPolicy::initSolver ( solverList );
    M_rhs.reset ( new vector_Type ( *M_solutionMap, Unique ) );
}

template< class AssemblyPolicy, class SolverPolicy >
void
TimeIterationPolicyNonlinear<AssemblyPolicy, SolverPolicy>::
iterate ( vectorPtr_Type solution,
          bcContainerPtr_Type bchandler,
          const Real& currentTime )
{
    int subiter = 0;

    Real normRhs ( 0.0 );
    Real nonLinearResidual ( 0.0 );
    Real rhsIterNorm ( 0.0 );

    do
    {
        //
        // STEP 1: Updating the system
        //
        displayer().leaderPrint ( "Updating the system... " );
        *M_rhs = 0.0;
        M_systemMatrix.reset ( new matrix_Type ( *M_solutionMap ) );
        AssemblyPolicy::assembleSystem ( M_systemMatrix, M_rhs, solution, SolverPolicy::preconditioner() );
        displayer().leaderPrint ( "done\n" );

        //
        // STEP 2: Applying the boundary conditions
        //
        displayer().leaderPrint ( "Applying BC... " );
        bcManage ( *M_systemMatrix, *M_rhs, *uFESpace()->mesh(), uFESpace()->dof(), *bchandler, uFESpace()->feBd(), 1.0, currentTime );
        M_systemMatrix->globalAssemble();
        displayer().leaderPrint ( "done\n" );

        // Norm of the rhs needed for the nonlinear convergence test
        if ( subiter == 0 )
        {
            normRhs = M_rhs->norm2();
        }

        //
        // STEP 3: Computing the residual
        //

        // Computing the RHS as RHS=b-Ax_k
        vector_Type Ax ( solution->map() );
        M_systemMatrix->matrixPtr()->Apply ( solution->epetraVector(), Ax.epetraVector() );

        Ax.epetraVector().Update (-1, M_rhs->epetraVector(), 1);
        nonLinearResidual = Ax.norm2();

        displayer().leaderPrint ( "Nonlinear residual          : ", nonLinearResidual, "\n" );
        displayer().leaderPrint ( "Nonlinear residual (scaled) : ", nonLinearResidual / normRhs, "\n" );

        if ( nonLinearResidual > M_nonLinearTolerance * normRhs )
        {
            displayer().leaderPrint ( "---\nSubiteration [", ++subiter, "]\n" );

            // Extra information if we want to know the exact residual
            if ( M_computeResidual )
            {
                rhsIterNorm = M_rhs->norm2();
            }

            //
            // Solving the system
            //
            displayer().leaderPrint ( "Solving the system... \n" );
            *solution = 0.0;
            SolverPolicy::solve ( M_systemMatrix, M_rhs, solution );
            // int numIter = SolverPolicy::solve( M_systemMatrix, M_rhs, solution );
            // numIterSum += numIter; //

            if ( M_computeResidual )
            {
                vector_Type Ax ( solution->map() );
                vector_Type res ( *M_rhs );
                M_systemMatrix->matrixPtr()->Apply ( solution->epetraVector(), Ax.epetraVector() );
                res.epetraVector().Update ( -1, Ax.epetraVector(), 1 );
                Real residual;
                res.norm2 ( &residual );
                residual /= rhsIterNorm;
                displayer().leaderPrint ( "Scaled residual: ", residual, "\n" );
            }
        }
    }
    while ( nonLinearResidual > M_nonLinearTolerance * normRhs );

    displayer().leaderPrint ( "Nonlinear iterations           : ", subiter, "\n" );
}

} // namespace LifeV

#endif /* TIMEITERATIONPOLICYNONLINEAR_HPP */
