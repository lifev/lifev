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
   \file SolverUMFPACK.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-04
 */
#include <debug.hpp>

#include <SolverUMFPACK.hpp>



namespace LifeV
{
void
SolverUMFPACK::reportInfo()
{
    umfpack_dl_report_info( _M_Control.get(), _M_Info.get() );
}
void
SolverUMFPACK::reportStatus( int status )
{
    umfpack_dl_report_status( _M_Control.get(), status);
}
void SolverUMFPACK::prepareSolve()
{
    if ( _M_matrix_reset )
    {
        if ( _M_symbolic )
        {
            Debug(5100) << "[SolverUMFPACK::prepareSolve] Destroying symbolic factorization\n";

            umfpack_dl_free_symbolic( &_M_symbolic );
            _M_symbolic = 0;
        }
        Debug(5100) << "[SolverUMFPACK::prepareSolve] computing symbolic factorization\n";
        int status = umfpack_dl_symbolic( _M_nrows,
                                          _M_nrows,
                                          _M_ia,
                                          _M_ja,
                                          _M_v,
                                          &_M_symbolic,
                                          _M_Control.get(),
                                          _M_Info.get() );
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
            umfpack_dl_free_numeric( &_M_numeric );
            _M_numeric = 0;
        }
        Debug(5100) << "[SolverUMFPACK::prepareSolve] computing numeric factorization\n";
        int status = umfpack_dl_numeric( _M_ia,
                                         _M_ja,
                                         _M_v,
                                         _M_symbolic, &_M_numeric,
                                         _M_Control.get(),
                                         _M_Info.get() );
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
