/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-08-29

  Copyright (C) 2004 EPFL

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
   \file SolverPETSC.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-08-29
*/
#ifndef __SolverPETSC_H
#define __SolverPETSC_H 1

#include <iostream>

extern "C"
{
#if defined(HAVE_PETSC_H)
#include <petsc.h>
#include <petscksp.h>
#include <petscpc.h>
#endif /* HAVE_PETSC_H */
}

#include <vecUnknown.hpp>

namespace LifeV
{
/*!
  \typedef enum SMonitorType
  Define the monitor type for the Krylov Space Solvers
*/
typedef enum PetscMonitorType
{
    KSP_NO_MONITOR = 0,  /**< no monitor */
    KSP_DEFAULT_MONITOR,  /**< preconditionned error monitor */
    KSP_TRUE_MONITOR     /**< preconditionned error and true error monitor */
};

/*!
  \class SolverPETSC
  \brief wrap petsc linear solvers

  You should have a look at PETSC documentation for further details.

  By default the solver is gmres and the preconditioner is ilu.

  Here is the list of solvers available:

  -# "richardson"
  -# "chebychev"
  -# "cg"
  -# "gmres"
  -# "tcqmr"
  -# "bcgs"
  -# "cgs"
  -# "tfqmr"
  -# "cr"
  -# "lsqr"
  -# "preonly"
  -# "qcg"
  -# "bicg"
  -# "fgmres"
  -# "minres"
  -# "symmlq"
  -# "lgmres"

  Here is the list of preconditioners available:

  -# "none"
  -# "jacobi"
  -# "sor"
  -# "lu"
  -# "shell"
  -# "bjacobi"
  -# "mg"
  -# "eisenstat"
  -# "ilu"
  -# "icc"
  -# "asm"
  -# "sles"
  -# "composite"
  -# "redundant"
  -# "spai"
  -# "milu"
  -# "nn"
  -# "cholesky"
  -# "ramg"
  -# "samg"
  -# "pbjacobi"
  -# "multilevel"
  -# "schur"
  -# "esi"
  -# "petscesi"
  -# "mat"
  -# "hypre"

  @author Christophe Prud'homme
  @see
*/
class SolverPETSC
{
public:


    /** @name Typedefs
     */
    //@{

    typedef double value_type;
    typedef SolverPETSC solver_type;
    //typedef St::Array::SArray<double,1>::type array_type;
    typedef Vector array_type;

    class Private;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    /*!
      The solver and preconditionner are the ones defined in petsc/petscksp.h
      \arg ksp krylov subspace method
      \arg pc preconditionner
    */
    SolverPETSC( std::string const& ksp = "gmres",
                 std::string const& pc = "ilu",
                 PetscMonitorType monitor = KSP_DEFAULT_MONITOR );

    //! create a new instance
    static SolverPETSC* New();

    //! destructor
    ~SolverPETSC();

    //@}

    /** @name Accessors
     */
    //@{

    double residualNorm() const;

    //! get the petsc preconditioner
    PC const& preconditioner() const;

    //! get the petsc Krylov solver
    KSP const& krylovSolver() const;

    //@}

    /** @name  Mutators
     */
    //@{

    void setMatrix( uint, const uint*, const uint*, const double* );
    void setMatrixTranspose( uint, const uint*, const uint*, const double* );
    void setTolerances( double = PETSC_DEFAULT, double = PETSC_DEFAULT, double = PETSC_DEFAULT, int = PETSC_DEFAULT );

    //! get the petsc preconditioner
    PC & preconditioner();

    //! get the petsc Krylov solver
    KSP & krylovSolver();

    //@}

    /** @name  Methods
     */
    //@{

    /*
      solve the problem \f$ A x = b \f$

      \c A has been entered via \c setMatrix .

      \c __ptype can have the following values :

      -# \c SAME_PRECONDITIONER -
      Pmat is identical during successive linear solves.
      This option is intended for folks who are using
      different Amat and Pmat matrices and want to reuse the
      same preconditioner matrix.  For example, this option
      saves work by not recomputing incomplete factorization
      for ILU/ICC preconditioners.


      -# \c SAME_NONZERO_PATTERN :
      Pmat has the same nonzero structure during
      successive linear solves.

      -# \c DIFFERENT_NONZERO_PATTERN -
      Pmat does not have the same nonzero structure.
    */
    int solve( array_type& __X, array_type const& __B, MatStructure __ptype = SAME_NONZERO_PATTERN );

    /*
      solve the transpose problem \f$ A^T x = b  \f$

      \c A has been entered via \c setMatrix .

      \c __ptype can have the following values :

      -# \c SAME_PRECONDITIONER -
      Pmat is identical during successive linear solves.
      This option is intended for folks who are using
      different Amat and Pmat matrices and want to reuse the
      same preconditioner matrix.  For example, this option
      saves work by not recomputing incomplete factorization
      for ILU/ICC preconditioners.


      -# \c SAME_NONZERO_PATTERN :
      Pmat has the same nonzero structure during
      successive linear solves.

      -# \c DIFFERENT_NONZERO_PATTERN -
      Pmat does not have the same nonzero structure.
    */
    int solveTranspose( array_type& __X, array_type const& __B, MatStructure __ptype = SAME_NONZERO_PATTERN );


    //@}

private:
    Private* _M_p;

};
}
#endif /* __SolverPETSC_H */
