/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-04

  Copyright (C) 2004 EPFL

  This library is free software; you can redlstribute it and/or
  modlfy it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is dlstributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file SolverUMFPACK.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-04
 */
#include <life/lifecore/debug.hpp>

#include <life/lifealg/SolverUMFPACK.hpp>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace LifeV
{
class SolverUMFPACK::Pimpl
{
public:

    typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::column_major> matrix_type;

    Pimpl()
        :
        _M_mat()
        {}

    Pimpl( Pimpl const& __p )
        :
        _M_mat( __p._M_mat )
        {}


    matrix_type _M_mat;
};
SolverUMFPACK::SolverUMFPACK()
    :
    _M_p( new Pimpl ),
    _M_matrix_reset( true ),
    _M_matrix_values_reset( true ),
    _M_symbolic( 0 ),
    _M_numeric( 0 ),
    _M_Control( new double[UMFPACK_CONTROL] ),
    _M_Info( new double[UMFPACK_INFO] )
{
    ::umfpack_di_defaults( _M_Control );
    _M_Control[UMFPACK_PRL] = 6;
    //_M_Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
    //_M_Control[UMFPACK_AMD_DENSE] = 1;
}
SolverUMFPACK::SolverUMFPACK( SolverUMFPACK const & umfpackSolver )
    :
    _M_p( umfpackSolver._M_p ),
    _M_matrix_reset( umfpackSolver._M_matrix_reset ),
    _M_matrix_values_reset( umfpackSolver._M_matrix_values_reset ),
    _M_symbolic( 0 ),
    _M_numeric( 0 ),
    _M_Control( umfpackSolver._M_Control ),
    _M_Info( umfpackSolver._M_Info )
{

}
SolverUMFPACK::~SolverUMFPACK()
{
    if ( _M_numeric )
        ::umfpack_di_free_numeric( &_M_numeric );

    if ( _M_symbolic )
        ::umfpack_di_free_symbolic( &_M_symbolic );

}

void
SolverUMFPACK::reportInfo()
{
    umfpack_di_report_info( _M_Control, _M_Info );
}
void
SolverUMFPACK::reportStatus( int status )
{
    umfpack_di_report_status( _M_Control, status);
}
void
SolverUMFPACK::setMatrix( const matrix_type& m,
                          bool samePattern )
{
    _M_p->_M_mat = m;
    _M_matrix_values_reset = true;
    if ( !samePattern )
        _M_matrix_reset = true;
}
void
SolverUMFPACK::setMatrix( const CSRMatr<CSRPatt, value_type>& m,
                          bool samePattern )
{
    Debug( 5100 ) << "copying csr matrix into a csc matrix\n";
    Debug( 5100 ) << "nrows = " << m.Patt()->nRows() << "\n";
    Debug( 5100 ) << "ncols = " << m.Patt()->nCols() << "\n";
    Debug( 5100 ) << "nnz   = " << m.Patt()->nNz() << "\n";
    _M_p->_M_mat.resize( m.Patt()->nRows(), m.Patt()->nRows(), false );

    CSRPatt::Container::const_iterator ia = m.Patt()->give_ia().begin();

    UInt nrows = m.Patt()->nRows();
    for ( UInt iRow = 0; iRow < nrows; ++iRow, ++ia )
    {
        for ( UInt i = *ia - OFFSET; i < *( ia + 1 ) - OFFSET; ++i )
        {
            UInt iCol = m.Patt()->give_ja()[ i ] - OFFSET;
            _M_p->_M_mat( iRow,  iCol ) = m.get_value( iRow, iCol );
        }
    }
    _M_matrix_values_reset = true;
    if ( !samePattern )
        _M_matrix_reset = true;
}
void
SolverUMFPACK::solve( array_type& __X, array_type const& __B )
{

    prepareSolve();

    Debug(5100) << "[SSolverUMFPACK:;solve] solve A x = b using UMFPACK version " << UMFPACK_VERSION << "\n";
    int status = umfpack_di_solve( UMFPACK_A,
                                   ( const int* )&_M_p->_M_mat.index1_data()[0],
                                   ( const int* )&_M_p->_M_mat.index2_data()[0],
                                   &_M_p->_M_mat.value_data()[0],
                                   boost::addressof( __X[0] ),
                                   boost::addressof( __B[0] ),
                                   _M_numeric,
                                   _M_Control,
                                   _M_Info );
    if (status != UMFPACK_OK)
    {
        reportInfo();
        reportStatus( status );
        Error() << "[SSolverUMFPACK::solve] solve failed\n";
    }
}
void SolverUMFPACK::prepareSolve()
{
    if ( _M_matrix_reset )
    {
        if ( _M_symbolic )
        {
            Debug(5100) << "[SolverUMFPACK::prepareSolve] Destroying symbolic factorization\n";

            umfpack_di_free_symbolic( &_M_symbolic );
            _M_symbolic = 0;
        }

        Debug(5100) << "[SolverUMFPACK::prepareSolve] computing symbolic factorization\n";
        int status = umfpack_di_symbolic( _M_p->_M_mat.size1(),
                                          _M_p->_M_mat.size2(),
                                          ( const int* )&_M_p->_M_mat.index1_data()[0],
                                          ( const int* )&_M_p->_M_mat.index2_data()[0],
                                          &_M_p->_M_mat.value_data()[0],
                                          &_M_symbolic,
                                          _M_Control,
                                          _M_Info );
        if (status != UMFPACK_OK)
        {
            reportInfo();
            reportStatus( status );
            Error() << "[SolverUMFPACK::prepareSolve] symbolic factorization failed\n";
        }
    }
    if ( _M_matrix_reset || _M_matrix_values_reset )
    {
        if ( _M_numeric )
        {
            Debug(5100) << "[SolverUMFPACK::prepareSolve] Destroying numeric factorization\n";
            umfpack_di_free_numeric( &_M_numeric );
            _M_numeric = 0;
        }
        Debug(5100) << "[SolverUMFPACK::prepareSolve] computing numeric factorization\n";
        int status = umfpack_di_numeric( ( const int* )&_M_p->_M_mat.index1_data()[0],
                                         ( const int* )&_M_p->_M_mat.index2_data()[0],
                                         &_M_p->_M_mat.value_data()[0],
                                         _M_symbolic, &_M_numeric,
                                         _M_Control,
                                         _M_Info );
        if (status != UMFPACK_OK)
        {
            reportInfo();
            reportStatus( status );
            Error() << "[SolverUMFPACK::prepareSolve] numeric factorization failed\n";
        }
    }
    _M_matrix_reset = false;
    _M_matrix_values_reset = false;
}

}
