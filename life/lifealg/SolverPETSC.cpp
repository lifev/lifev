/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-
 
  This file is part of the LifeV libraries.
 
  Author: Christophe Prud'homme <christophe.prudhomme@epfl.ch>
 
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
/** \file SolverPETSC.cpp
    Implementation file for the LifeV-PETSC interface.
    
*/
#if defined(__GNUG__)
#pragma implementation
#endif /* __GNUG__ */

#include <SolverPETSC.hpp>

extern "C"
{
#include <petsc.h>
#include <petscmat.h>
#include <petscsles.h>
};


namespace life
{

    class SolverPETSC::Private
    {
    public:

        Private()
            :
            __tolerance( 1e-11 ),
            _M_use_A( false ),
            _M_use_A_t( false )
            {}
        double __tolerance;
        Mat __A;
        Mat __A_t;

        SLES __sles;
        PC __pc;
        KSP __ksp;

        mutable bool _M_use_A;
        mutable bool _M_use_A_t;

    };



    SolverPETSC::SolverPETSC( std::string const& __ksp_type,
                              std::string const& __pc_type,
                              SMonitorType __monitor )
        :
        _M_p ( new Private )
    {
        int ierr = SLESCreate( PETSC_COMM_WORLD, &_M_p->__sles ); //CHKERRQ(ierr);

        ierr = SLESGetKSP( _M_p->__sles, & _M_p->__ksp ); //CHKERRQ(ierr);

        ierr = SLESGetPC( _M_p->__sles, &_M_p->__pc ); //CHKERRQ(ierr);

        ierr = KSPSetType( _M_p->__ksp, const_cast<char*> ( __ksp_type.c_str() ) ); //CHKERRQ(ierr);

        //! set the monitor type
        switch ( __monitor )
        {
            case KSP_DEFAULT_MONITOR:
                ierr = KSPSetMonitor( _M_p->__ksp, &KSPDefaultMonitor, 0, 0 );
                break;

            case KSP_TRUE_MONITOR:
                ierr = KSPSetMonitor( _M_p->__ksp, &KSPTrueMonitor, 0, 0 );
                break;

            case KSP_NO_MONITOR:
            default:
                // do nothing
                break;
        }
        ierr = PCSetType( _M_p->__pc, const_cast<char*> ( __pc_type.c_str()  ) ); //CHKERRQ(ierr);

        if ( __pc_type == "ilu" )
        {
            ierr = PCILUSetLevels( _M_p->__pc, 2 ); //CHKERRQ(ierr);
        }
        if ( __pc_type == "icc" )
        {
            ierr = PCICCSetLevels( _M_p->__pc, 2 ); //CHKERRQ(ierr);
        }
        ierr = KSPSetTolerances( _M_p->__ksp, _M_p->__tolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT ); //CHKERRQ(ierr);

        ierr = SLESSetFromOptions( _M_p->__sles ); //CHKERRQ(ierr);
    }

    SolverPETSC::~SolverPETSC()
    {
        int __ierr;
        if ( _M_p->_M_use_A )
        {
            __ierr = MatDestroy( _M_p->__A ); ////CHKERRQ( __ierr );
        }
        if ( _M_p->_M_use_A_t )
        {
            __ierr = MatDestroy( _M_p->__A_t ); ////CHKERRQ( __ierr );
        }
    }

    SolverPETSC*
    SolverPETSC::New()
    {
        return new SolverPETSC;
    }

    /*!
      \brief Gets the last (approximate preconditioned) residual norm that has been computed. 
    
      \return last (approximate preconditioned) residual norm
    */
    double
    SolverPETSC::residualNorm() const
    {
        double __residual;
        KSPGetResidualNorm( _M_p->__ksp, &__residual );
        return __residual;
    }
    
    PC const& SolverPETSC::preconditioner() const { return _M_p->__pc; }
    PC &      SolverPETSC::preconditioner()       { return _M_p->__pc; }
    
    KSP const& SolverPETSC::krylovSolver() const { return _M_p->__ksp; }
    KSP &      SolverPETSC::krylovSolver()       { return _M_p->__ksp; }
    
    
    /*!
      \brief Sets the relative, absolute, divergence, and maximum iteration tolerances 
      used by the default KSP convergence testers. 
    
      Use PETSC_DEFAULT to retain the default value of any of the tolerances. 
    
      See KSPDefaultConverged() for details on the use of these parameters in the default convergence test. 
      See also KSPSetConvergenceTest() for setting user-defined stopping criteria. 
    
      \arg rtol - the relative convergence tolerance (relative decrease in the residual norm) 
      \arg atol - the absolute convergence tolerance (absolute size of the residual norm) 
      \arg dtol - the divergence tolerance (amount residual can increase before KSPDefaultConverged() concludes that the method is diverging) 
      \arg maxits - maximum number of iterations to use 
    */
    void
    SolverPETSC::setTolerances( double __rtol, double __atol, double __dtol, int __maxits  )
    {
        KSPSetTolerances( _M_p->__ksp, __rtol, __atol, __dtol, __maxits  );
    }

    void
    SolverPETSC::setMatrix( uint __nrows, uint* __r, uint *__i, double* __v )
    {
        MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,
                                  __nrows, __nrows,
                                  (int*)__r,
                                  (int*)__i,
                                  __v,
                                  &_M_p->__A );

        _M_p->_M_use_A = true;
    }

    void
    SolverPETSC::setMatrixTranspose( uint __nrows, uint* __r, uint *__i, double* __v )
    {
        MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,
                                  __nrows, __nrows,
                                  (int*)__r,
                                  (int*)__i,
                                  __v,
                                  //const_cast<int*>( __M.getRowIndices() ),
                                  //const_cast<int*>( __M.getMatrixIndices() ),
                                  //const_cast<double*>( __M.getMatrixValues() ),
                                  &_M_p->__A_t );
        _M_p->_M_use_A_t = true;

    }

    //! solve \f[ A X = B \f]
    /*!
    
    \param __X  the solution
    \param __B	the right hand side
    \return the number of iterations
    */
    int
    SolverPETSC::solve( array_type& __X, array_type const& __B, MatStructure  __ptype )
    {
        int __rowA;
        int __colA;
        MatGetSize(_M_p->__A,&__rowA,&__colA); //CHKERRQ(__ierr);

        int __iterations;

        std::cerr << "[SolverPETSC::solve]  Solving primal\n";
        Vec __x;
        VecCreateSeqWithArray( PETSC_COMM_SELF, __X.size(), &__X[0], &__x );

        Vec __b;
        VecCreateSeqWithArray( PETSC_COMM_SELF, __B.size(), &__B[0], &__b );

    
        SLESSetOperators( _M_p->__sles, _M_p->__A, _M_p->__A, __ptype  ); //CHKERRQ(__ierr);
        SLESSolve( _M_p->__sles, __b, __x , &__iterations ); //CHKERRQ(__ierr);

        /*
          View solver info; we could instead use the option -sles_view to
          print this info to the screen at the conclusion of SLESSolve().
        */
        SLESView( _M_p->__sles, PETSC_VIEWER_STDOUT_WORLD );

        // Petsc won't deallocate the memory so __X and __B contains the informations
        VecDestroy( __x );
        VecDestroy( __b );

        std::cerr << "[SolverPETSC::solve] Solving primal done in " << __iterations << " iterations\n";
        return __iterations;
    }

    int
    SolverPETSC::solveTranspose( array_type& __X, array_type const& __B, MatStructure __ptype )
    {
        int __rowA;
        int __colA;
        MatGetSize(_M_p->__A_t,&__rowA,&__colA); ////CHKERRQ(__ierr);

        int __iterations;

        std::cerr << "[SolverPETSC::solveTranspose] Solving transpose\n";
        Vec __x;
        VecCreateSeqWithArray( PETSC_COMM_SELF, __X.size(), &__X[0], &__x );

        Vec __b;
        VecCreateSeqWithArray( PETSC_COMM_SELF, __B.size(), &__B[0], &__b );

        SLESSetOperators( _M_p->__sles, _M_p->__A_t, _M_p->__A_t, __ptype ); ////CHKERRQ(__ierr);
        SLESSolve( _M_p->__sles, __b, __x , &__iterations ); ////CHKERRQ(__ierr);

        // Petsc won't deallocate the memory so __X and __B contains the informations
        VecDestroy( __x );
        VecDestroy( __b );

        std::cerr << "[SolverPETSC::solveTranspose] Solving transpose done in " << __iterations << " iterations\n";
        return __iterations;
    }

}
