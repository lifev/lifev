/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
            Christoph Winkelmann <christoph.winkelmann@epfl.ch>
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
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-08-29
*/
#ifndef __SolverPETSC_H
#define __SolverPETSC_H 1

#include <iostream>

#if defined(HAVE_PETSC_H)
#define PETSC_USE_EXTERN_CXX
#include <petsc.h>
#include <petscksp.h>
#include <petscpc.h>
#endif /* HAVE_PETSC_H */

#include <singleton.hpp>

#include <vecUnknown.hpp>
#include <sparseArray.hpp>

class GetPot;

namespace LifeV
{
/*!
  \typedef enum SMonitorType
  Define the monitor type for the Krylov Space Solvers
*/
typedef enum PetscMonitorType
{
    KSP_NO_MONITOR = 0,   /**< no monitor */
    KSP_DEFAULT_MONITOR,   /**< preconditionned error monitor */
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

  \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
  \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
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
      \param ksp krylov subspace method
      \param pc preconditionner
      \param monitor Petsc monitor called at each iteration
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

    /*!
      \brief Gets the last (approximate preconditioned) residual norm that has
      been computed.

      \return last (approximate preconditioned) residual norm
    */
    double residualNorm() const;

    /*!
      \brief Gets a condition number estimate of the preconditioned system.

      The estimate is the quotient of the minimal and the maximal singular
      value found in the Krylov space used in the last solve.

      Note: The necessary data is only computed during the solve if the option
      ksp_set_compute_singular_values is set to true in the data file.

      \return estimated condition number
    */
    double condEst() const;

    //! Gets the number of iterations performed in the last solve.
    int iterations() const;

    //! Returns whether the last solve converged
    bool converged() const;

    //! get the petsc preconditioner
    PC const& preconditioner() const;

    //! get the petsc Krylov solver
    KSP const& krylovSolver() const;

    //@}

    /** @name  Mutators
     */
    //@{

    //! set matrix from raw CSR arrays
    void setMatrix( uint, const uint*, const uint*, const double* );

    //! set matrix from CSRMatr
    template <typename PatternType>
    void setMatrix( const CSRMatr<PatternType, value_type>& m )
    {
        _M_tempPattern.reset( 0 );
        _M_tempMatrix.reset( 0 );
        setMatrix( m.Patt() ->nRows(),
                   m.Patt() ->giveRawCSR_ia(),
                   m.Patt() ->giveRawCSR_ja(),
                   m.giveRawCSR_value() );
    }

    /** set matrix from MSRMatr
     *
     *  Warning: The matrix is converted to CSR. This method provides ease of
     *  use, possibly for the sake of efficiency.
     */
    void setMatrix( const MSRMatr<value_type>& m );

    void setMatrixTranspose( uint, const uint*, const uint*, const double* );
    /*!
      \brief Sets the relative, absolute, divergence, and maximum iteration
      tolerances used by the default KSP convergence testers.

      Use PETSC_DEFAULT to retain the default value of any of the tolerances.

      See KSPDefaultConverged() for details on the use of these parameters in
      the default convergence test. See also KSPSetConvergenceTest() for
      setting user-defined stopping criteria.

      \arg rtol - the relative convergence tolerance (relative decrease in the
      residual norm)
      \arg atol - the absolute convergence tolerance (absolute size of the
      residual norm)
      \arg dtol - the divergence tolerance (amount residual can increase before
      KSPDefaultConverged() concludes that the method is diverging)
      \arg maxits - maximum number of iterations to use
    */
    void setTolerances( double = PETSC_DEFAULT,
                        double = PETSC_DEFAULT,
                        double = PETSC_DEFAULT,
                        int = PETSC_DEFAULT );

    /*!
      \brief set the null space of the matrix to be solved for
      
      \param nullSpace orthonormal basis of the null space
    */
    void setNullSpace( const std::vector<const Vector*>& nullSpace );
    
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
    void solve( array_type& __X,
                array_type const& __B,
                MatStructure __ptype = SAME_PRECONDITIONER );

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
    void solveTranspose( array_type& __X,
                         array_type const& __B,
                         MatStructure __ptype = SAME_PRECONDITIONER );


    //@}

    /*!
      @brief Adds options from data file to PETSC database and sets these new
      options for this solver.

      Any option of the form NAME = VALUE is passed to PETSC as command line
      option -NAME VALUE.

      Example: ksp_type = gmres. See the PETSC documentation for more
      available options.

      In addition, there are the following options:

      @arg nokspview = false | true - Do not show parameters of Krylov space
      solver. Default: false
      @arg quiet = false | true - Do not show log messages. Only error
      messages are printed to standard error. Default: false
      @arg ksp_set_compute_singular_values = false | true - Calculate singular
      values when solving. Needed if condEst shall not return NaN. Default:
      false

      @param dataFile GetPot object containing the options from the data file
      @param section section in the GetPot object containing the PETSC stuff

      You should have a look at PETSC documentation for further details.
      @author Christoph Winkelmann
      @see http://www.mcs.anl.gov/petsc/
     */
    void setOptionsFromGetPot( const GetPot& dataFile,
                               std::string section = "petsc" );

private:
    std::auto_ptr<Private> _M_p;

    //! CSRPatt converted from MSRPatt if matrix given as MSRMatr
    std::auto_ptr<CSRPatt> _M_tempPattern;

    //! CSRMatr converted from MSRMatr if given as such
    std::auto_ptr<CSRMatr<CSRPatt, value_type> > _M_tempMatrix;

    //! common code for solve and solveTranspose
    void _F_solveCommon( Mat const& __A,
                         array_type& __X,
                         array_type const& __B,
                         MatStructure __ptype = SAME_NONZERO_PATTERN);
};

/*!
  \class PETSCManager
  \brief initialization and finalization for PETSC linear solvers

  Should be used as PETSC, which is a typedef to LifeV::singleton<PETSCManager>

  @author Christoph Winkelmann
  @see http://www.mcs.anl.gov/petsc/
*/
class PETSCManager
{
public:

    PETSCManager()
        :
        _M_initialized( false )
        {
        }
    /**
       initialized PETSC
    */
    bool initialize()
        {
            if ( !isInitialized() )
            {
                Debug( 7010 ) << "initialize PETSC with dummy arguments\n";

                // dummy arguments for PetscInitialize
                int __argc = 1;
                char** __argv = ( char** )malloc( 2*sizeof( char* ) );
                __argv[0] = ( char* )malloc( 2*sizeof( char ) );
                __argv[0][0] = 't';
                __argv[0][1] = '\0';

                // needed to avoid segmentation fault in PETSc 2.2.1
                __argv[1] = 0;

                _M_initialized = initialize( __argc, __argv );

                free( __argv[0] );
                free( __argv );

            }
            return _M_initialized;
        }

    /**
       initialized PETSC
    */
    bool initialize( int __argc, char** __argv )
        {
            if ( !isInitialized() )
            {
                Debug( 7010 ) << "initialize PETSC\n";
                //MPI_Init( &__argc, &__argv );
                PetscInitialize( &__argc, &__argv,PETSC_NULL,PETSC_NULL);

                _M_initialized = true;

            }
            return _M_initialized;
        }

    //! return \c true if PETSC manager is initialized, \c false otherwise
    bool isInitialized() const
        {
            return _M_initialized;
        }

    //! finalize PETSC environment
    void finalize()
        {
            if ( isInitialized() )
            {
                //MPI_Finalize();
                PetscFinalize();
                _M_initialized = false;
            }
        }
    //! Finalizes PETSC
    ~PETSCManager()
        {
            if ( isInitialized() )
                finalize();
        }
private:
    bool _M_initialized;
};
typedef singleton<PETSCManager> PETSC;
}
#endif /* __SolverPETSC_H */
