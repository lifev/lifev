/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-04

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
   \file SolverUMFPACK.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-04
 */
#ifndef __SolverUMFPACK_H
#define __SolverUMFPACK_H 1

extern "C"
{
#include <umfpack.h>
};

#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <vecUnknown.hpp>

namespace LifeV
{
/*!
  \class SolverUMFPACK
  \brief Interface for the UMFPACK Solver

  UMFPACK is a direct Solver for (un)symmetric problem \f$ A x = b \f$.

  @author Christophe Prud'homme
*/
class SolverUMFPACK
{
public:


    /** @name Typedefs
     */
    //@{

    typedef double value_type;
    typedef Vector array_type;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    SolverUMFPACK()
        :
        _M_nrows( 0 ),
        _M_ia( 0 ),
        _M_ja( 0 ),
        _M_v( 0 ),
        _M_matrix_reset( true ),
        _M_matrix_values_reset( true ),
        _M_symbolic( 0 ),
        _M_numeric( 0 ),
        _M_Control( new double[UMFPACK_CONTROL] ),
        _M_Info( new double[UMFPACK_INFO] )
        {
            ::umfpack_dl_defaults( _M_Control.get() );
        }
    SolverUMFPACK( SolverUMFPACK const & umfpackSolver )
        :
        _M_nrows( umfpackSolver._M_nrows ),
        _M_ia( umfpackSolver._M_ia ),
        _M_ja( umfpackSolver._M_ja ),
        _M_v( umfpackSolver._M_v ),
        _M_matrix_reset( umfpackSolver._M_matrix_reset ),
        _M_matrix_values_reset( umfpackSolver._M_matrix_values_reset ),
        _M_symbolic( 0 ),
        _M_numeric( 0 ),
        _M_Control( umfpackSolver._M_Control ),
        _M_Info( umfpackSolver._M_Info )
        {

        }
    ~SolverUMFPACK()
        {
            if ( _M_numeric )
                ::umfpack_dl_free_numeric( &_M_numeric );

            if ( _M_symbolic )
                ::umfpack_dl_free_symbolic( &_M_symbolic );

            delete[] _M_Control;
            delete[] _M_Info;
        }

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


    //@}

    /** @name  Mutators
     */
    //@{

    //! set matrix from raw CSR arrays
    void setMatrix( uint __N, const uint* __ia, const uint* __ja, const double* __v )
        {
            _M_nrows = __N;
            _M_ia = ( long int* )__ia;
            _M_ja = ( long int* )__ja;
            _M_v = __v;
        }

    //@}

    /** @name  Methods
     */
    //@{


    //! solve A X = B
    /*!

    \param __X  the solution
    \param __B	the right hand side
    \return the number of iterations
    */
    void solve( array_type& __X, array_type const& __B )
        {

            prepareSolve();

            Debug(5100) << "[SSolverUMFPACK:;solve] solve A x = b using UMFPACK version " << UMFPACK_VERSION << "\n";
            int status = umfpack_dl_solve( UMFPACK_A,
                                           _M_ia,
                                           _M_ja,
                                           _M_v,
                                           boost::addressof( __X[0] ),
                                           boost::addressof( __B[0] ),
                                           _M_numeric,
                                           _M_Control.get(),
                                           _M_Info.get() );

            if (status != UMFPACK_OK)
            {
                reportInfo();
                reportStatus( status );
                Error() << "[SSolverUMFPACK::solve] solve failed\n";
            }
        }

    //! report some info about umfpack
    void reportInfo();


    /**
       report status of umfpack

       \param status status integer returned by umfpack routines
     */
    void reportStatus( int );

    //@}

private:

    void prepareSolve();

private:

    size_t _M_nrows;
    long int const* _M_ia;
    long int const* _M_ja;
    double const* _M_v;

    bool _M_matrix_reset;
    bool _M_matrix_values_reset;

    void *_M_symbolic;
    void *_M_numeric;

    boost::shared_ptr<double> _M_Control;
    boost::shared_ptr<double> _M_Info;

};
}
#endif /* __SolverUMFPACK_H */

