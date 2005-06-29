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
#include <lifeconfig.h>

#include <life/lifecore/debug.hpp>

#include <life/lifearray/vecUnknown.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/singleton.hpp>

#include <life/lifealg/SolverPETSC.hpp>

namespace LifeV
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

    PC __pc;
    KSP __ksp;

    mutable bool _M_use_A;
    mutable bool _M_use_A_t;

};



SolverPETSC::SolverPETSC( std::string const& __ksp_type,
                          std::string const& __pc_type,
                          PetscMonitorType __monitor )
        :
        _M_p ( new Private )
{
    int ierr = KSPCreate( PETSC_COMM_WORLD, &_M_p->__ksp );
    //CHKERRQ(ierr);

    ierr = KSPGetPC( _M_p->__ksp, &_M_p->__pc );
    //CHKERRQ(ierr);

    ierr = KSPSetType( _M_p->__ksp, const_cast<char*> ( __ksp_type.c_str() ) );
    //CHKERRQ(ierr);

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
    ierr = PCSetType( _M_p->__pc, const_cast<char*> ( __pc_type.c_str() ) );
    //CHKERRQ(ierr);

    if ( __pc_type == "ilu" )
    {
        //ierr = PCILUSetLevels( _M_p->__pc, 2 );
        //CHKERRQ(ierr);
        ierr = PCILUSetUseDropTolerance( _M_p->__pc, 1e-6, 0.1, 200 );
    }
    if ( __pc_type == "icc" )
    {
        ierr = PCICCSetLevels( _M_p->__pc, 2 );
        //CHKERRQ(ierr);
    }
    ierr = KSPSetTolerances( _M_p->__ksp, _M_p->__tolerance,
                             PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
    //CHKERRQ(ierr);

    ierr = KSPSetFromOptions( _M_p->__ksp );
    //CHKERRQ(ierr);
}

SolverPETSC::~SolverPETSC()
{
    int __ierr;
    if ( _M_p->_M_use_A )
    {
        __ierr = MatDestroy( _M_p->__A );
        //CHKERRQ( __ierr );
    }
    if ( _M_p->_M_use_A_t )
    {
        __ierr = MatDestroy( _M_p->__A_t );
        //CHKERRQ( __ierr );
    }
}

SolverPETSC*
SolverPETSC::New()
{
    return new SolverPETSC;
}

double
SolverPETSC::residualNorm() const
{
    double __residual;
    KSPGetResidualNorm( _M_p->__ksp, &__residual );
    return __residual;
}

double
SolverPETSC::condEst() const
{
    double __value;
    PetscTruth __computed;
    PetscTruth __found;
    PetscOptionsGetLogical( PETSC_NULL,
                            "-ksp_set_compute_singular_values",
                            &__computed, &__found );
    if ( true ) //__computed )
    {
        PetscReal __emin;
        PetscReal __emax;
        KSPComputeExtremeSingularValues( _M_p->__ksp, &__emax, &__emin);
        __value = __emax/__emin;
    }
    else
    {
        std::cerr << "Singular values have not been computed" << std::endl;
        std::cerr
            << "add 'ksp_set_compute_singular_values = true' to your data file"
            << std::endl;
        __value = std::numeric_limits<double>::signaling_NaN();
    }
    return __value;
}

int
SolverPETSC::iterations() const
{
    int __iter;
    KSPGetIterationNumber( _M_p->__ksp , &__iter );
    return __iter;
}

bool
SolverPETSC::converged() const
{
    KSPConvergedReason reason;
    KSPGetConvergedReason( _M_p->__ksp, &reason );
    return reason >= 0;
}

PC const& SolverPETSC::preconditioner() const
{
    return _M_p->__pc;
}
PC & SolverPETSC::preconditioner()
{
    return _M_p->__pc;
}

KSP const& SolverPETSC::krylovSolver() const
{
    return _M_p->__ksp;
}
KSP & SolverPETSC::krylovSolver()
{
    return _M_p->__ksp;
}

void
SolverPETSC::setTolerances( double __rtol, double __atol, double __dtol,
                            int __maxits )
{
    KSPSetTolerances( _M_p->__ksp, __rtol, __atol, __dtol, __maxits );
}

void
SolverPETSC::setMatrix( uint __nrows, const uint* __r, const uint *__i,
                        const double* __v )
{
    MatCreateSeqAIJWithArrays( PETSC_COMM_SELF,
                               __nrows, __nrows,
                               ( int* ) __r,
                               ( int* ) __i,
                               const_cast<PetscScalar*>( __v ),
                               &_M_p->__A );

    _M_p->_M_use_A = true;
}
void
SolverPETSC::setMatrix( size_t __nrows, const size_t* __r, const size_t *__i,
                        const double* __v )
{
    MatCreateSeqAIJWithArrays( PETSC_COMM_SELF,
                               __nrows, __nrows,
                               ( int* ) __r,
                               ( int* ) __i,
                               const_cast<PetscScalar*>( __v ),
                               &_M_p->__A );

    _M_p->_M_use_A = true;
}

void SolverPETSC::setMatrix( const MSRMatr<value_type>& m )
{
    _M_tempPattern.reset( new CSRPatt( *( m.Patt() ) ) );
    _M_tempMatrix.reset( new CSRMatr<CSRPatt, value_type>
                         ( *_M_tempPattern, m ));
    setMatrix( _M_tempPattern->nRows(),
               _M_tempPattern->giveRawCSR_ia(),
               _M_tempPattern->giveRawCSR_ja(),
               _M_tempMatrix->giveRawCSR_value() );
}

void
SolverPETSC::setMatrixTranspose( uint __nrows, const uint* __r,
                                 const uint *__i, const double* __v )
{
    MatCreateSeqAIJWithArrays( PETSC_COMM_SELF,
                               __nrows, __nrows,
                               ( int* ) __r,
                               ( int* ) __i,
                               const_cast<PetscScalar*>( __v ),
                               //const_cast<int*>( __M.getRowIndices() ),
                               //const_cast<int*>( __M.getMatrixIndices() ),
                               //const_cast<double*>( __M.getMatrixValues() ),
                               &_M_p->__A_t );
    _M_p->_M_use_A_t = true;

}

void
SolverPETSC::setNullSpace( const std::vector<const Vector*>& __nullSpace )
{
    UInt dimNullSpace = __nullSpace.size();
    UInt dimVec = __nullSpace[ 0 ]->size();
    Vec* nullSpaceVecs = new Vec[dimNullSpace];
    for ( UInt i=0; i<dimNullSpace; ++i )
    {
        VecCreateSeqWithArray( PETSC_COMM_SELF,
                               dimVec,
                               __nullSpace[ i ]->data().begin(),
                               nullSpaceVecs+i );
    }
    MatNullSpace nullSpace;
    MatNullSpaceCreate( PETSC_COMM_WORLD,
                        PETSC_FALSE,
                        dimNullSpace,
                        nullSpaceVecs,
                        &nullSpace );
#if PETSC_VERSION == 220
    PCNullSpaceAttach( _M_p->__pc, nullSpace );
#else
    KSPSetNullSpace( _M_p->__ksp, nullSpace );
#endif
    delete[] nullSpaceVecs;
}

void
SolverPETSC::solve( array_type& __X,
                    array_type const& __B,
                    MatStructure __ptype )
{
    PetscTruth __quiet;
    PetscTruth __found;
    PetscOptionsGetLogical( PETSC_NULL, "-quiet", &__quiet, &__found );
    if ( !__quiet )
    {
        std::cout << "[SolverPETSC::solve] Solving primal\n";
    }
    _F_solveCommon( _M_p->__A, __X, __B, __ptype );
}

void
SolverPETSC::solveTranspose( array_type& __X,
                             array_type const& __B,
                             MatStructure __ptype )
{
    PetscTruth __quiet;
    PetscTruth __found;
    PetscOptionsGetLogical( PETSC_NULL, "-quiet", &__quiet, &__found );
    if ( !__quiet )
    {
        std::cout << "[SolverPETSC::solveTranspose] Solving transpose\n";
    }
    _F_solveCommon( _M_p->__A_t, __X, __B, __ptype );
}

void
SolverPETSC::_F_solveCommon( Mat const& __A,
                             array_type& __X,
                             array_type const& __B,
                             MatStructure __ptype)
{
    int __rowA;
    int __colA;
    MatGetSize( __A, &__rowA, &__colA );
    //CHKERRQ(__ierr);

    Vec __x;
    VecCreateSeqWithArray( PETSC_COMM_SELF, __X.size(), &__X[ 0 ], &__x );

    Vec __b;
    VecCreateSeqWithArray( PETSC_COMM_SELF, __B.size(), &__B[ 0 ], &__b );

    KSPSetOperators( _M_p->__ksp, __A, __A, __ptype );
    //CHKERRQ(__ierr);

    KSPSetInitialGuessNonzero( _M_p->__ksp, PETSC_TRUE );

#if PETSC_VERSION == 220

    KSPSetRhs( _M_p->__ksp, __b );
    KSPSetSolution( _M_p->__ksp, __x );
    KSPSolve( _M_p->__ksp );
    //CHKERRQ(__ierr);
#else

    PetscTruth __quiet;
    PetscTruth __found;
    PetscOptionsGetLogical( PETSC_NULL, "-quiet", &__quiet, &__found );

    PetscInt its = 0;

    KSPSolve( _M_p->__ksp, __b, __x );
    KSPConvergedReason reason;
    KSPGetConvergedReason( _M_p->__ksp, &reason );
    if ( reason == KSP_DIVERGED_INDEFINITE_PC )
    {
        std::cerr << "\nDivergence because of indefinite preconditioner;\n";
    }
    else if ( reason < 0 )
    {
        std::cerr << "\nOther kind of divergence: this should not happen.\n";
    }
    else if ( !__quiet )
    {
        KSPGetIterationNumber( _M_p->__ksp , &its );
        std::cout << "\nConvergence in " << ( int ) its << " iterations.\n";
    }
    if ( !__quiet )
    {
        std::cout << "\n";
    }
#endif
    /*
      View info about the solver
    */
    PetscTruth __flag;
    PetscTruth __flagFound;
    PetscOptionsGetLogical( PETSC_NULL, "-nokspview", &__flag, &__flagFound );
    if ( !__flag )
    {
        KSPView( _M_p->__ksp, PETSC_VIEWER_STDOUT_WORLD );
    }

    // Petsc won't deallocate the memory so __X and __B contains the
    // informations
    VecDestroy( __x );
    VecDestroy( __b );
}

void SolverPETSC::setOptionsFromGetPot( const GetPot& dataFile,
                                        std::string section )
{
    // get variable names
    const getpot::StringVector vars = dataFile.get_variable_names();

    // build up vector of options
    getpot::StringVector petscOpts;
    petscOpts.push_back( "LifeV::SolverPETSC" );
    for ( getpot::StringVector::const_iterator it = vars.begin();
            it != vars.end(); ++it )
    {
        std::string::size_type p = it->find( section + "/" );
        if ( p != std::string::npos )
        { // one might use == 0 to be more strict
            std::string name = it->substr( p + section.size() + 1 );
            std::string value = dataFile( it->c_str(), "" );
            petscOpts.push_back( "-" + name );
            petscOpts.push_back( value );
        }
    }

    // build up argv and argc
    char** argv = new char * [ petscOpts.size() + 1 ];
    int argc = 0;
    for ( getpot::StringVector::const_iterator it = petscOpts.begin();
            it != petscOpts.end(); ++it )
    {
        argv[ argc ] = new char[ it->size() + 1 ];
        strcpy( argv[ argc++ ], it->c_str() );
    }
    argv[ argc ] = 0;
    //     for(int i=0; i<argc; ++i) {
    //         std::cout << argv[i] << std::endl;
    //     }

    // set options
    PetscOptionsInsert( &argc, &argv, 0 );
    KSPSetFromOptions( _M_p->__ksp ); //CHKERRQ(ierr);
    for ( int i = 0; i < argc; ++i )
    {
        delete[] argv[ i ];
    }
    delete[] argv;

    // set compute singular values from options
    // (not supported directly by PETSc)
    PetscTruth __compute;
    PetscTruth __found;
    PetscOptionsGetLogical( PETSC_NULL,
                            "-ksp_set_compute_singular_values",
                            &__compute, &__found );
    KSPSetComputeSingularValues( _M_p->__ksp, __compute );
}

static const bool __petsc_initialize = PETSC::instance().initialize();

} // namespace LifeV
