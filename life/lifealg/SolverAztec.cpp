/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
      Date: 2004-09-22

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
   \file SolverAztec.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-09-22
*/
#include <iostream>
#include <fstream>

#include <boost/utility.hpp>

#include <SolverAztec.hpp>

class GetPot;

namespace LifeV
{

UInt SolverAztec::S_solverNumber = 100;

SolverAztec::SolverAztec(std::string filename)
        : M_matrix( 0 ), M_precond( 0 ), M_tempPattern( 0 ), M_tempMatrix( 0 )
{
    M_dataOrg[ AZ_N_internal ] = 0;
    M_dataOrg[ AZ_N_border ] = 0;
    M_dataOrg[ AZ_N_external ] = 0;
    M_dataOrg[ AZ_N_neigh ] = 0;
    M_dataOrg[ AZ_name ] = S_solverNumber++;
    AZ_set_proc_config( M_procConfig, AZ_NOT_MPI );

    // let dataAztec set the defaults
    GetPot dataFile(filename.c_str());
    DataAztec dataAztec( dataFile, "aztec" );
    dataAztec.aztecOptionsFromDataFile( M_options, M_params );

    // use ilu by default
    M_options[ AZ_precond ] = AZ_dom_decomp;
    M_options[ AZ_subdomain_solve ] = AZ_ilu;
}

SolverAztec::~SolverAztec()
{
    if ( M_matrix != 0 )
    {
        AZ_matrix_destroy( &M_matrix );
    }
    if ( M_precond != 0 )
    {
        AZ_precond_destroy( &M_precond );
    }
}

SolverAztec* SolverAztec::New()
{
    return new SolverAztec;
}

/*!
  \brief Gets the last residual norm that has been computed.

  \return last residual norm
*/
double SolverAztec::residualNorm() const
{
    return M_status[ AZ_r ];
}

int SolverAztec::iterations() const
{
    return (int)M_status[ AZ_its ];
}

bool SolverAztec::converged() const
{
    return M_status[ AZ_why ] == AZ_normal;
}

void SolverAztec::setMatrix( MSRMatr<value_type> const& newMatrix )
{
    M_tempPattern.reset( 0 );
    M_tempMatrix.reset( 0 );
    F_setMatrix( newMatrix );
}

void SolverAztec::setMatrix( CSRMatr<CSRPatt, value_type> const& newMatrix )
{
    M_tempPattern.reset( new MSRPatt( *( newMatrix.Patt() ) ) );
    M_tempMatrix.reset( new MSRMatr<value_type>( *M_tempPattern, newMatrix ) );
    F_setMatrix( *M_tempMatrix );
}

void SolverAztec::
setMatrixFree( int nEq, void* data,
               void ( *matvec ) ( double*, double*, AZ_MATRIX_STRUCT*, int* ) )
{
    M_dataOrg[ AZ_N_internal ] = nEq;
    if ( M_matrix != 0 )
    {
        AZ_matrix_destroy( &M_matrix );
    }
    if ( M_precond != 0 )
    {
        AZ_precond_destroy( &M_precond );
    }
    M_matrix = AZ_matrix_create( nEq );
    AZ_set_MATFREE( M_matrix, data, matvec );
}

void SolverAztec::F_setMatrix( MSRMatr<value_type> const& newMatrix )
{
    int nEq = newMatrix.Patt() ->nRows();
    M_dataOrg[ AZ_N_internal ] = nEq;
    if ( M_matrix != 0 )
    {
        AZ_matrix_destroy( &M_matrix );
    }
    if ( M_precond != 0 )
    {
        AZ_precond_destroy( &M_precond );
    }
    M_matrix = AZ_matrix_create( nEq );
    M_precond = AZ_precond_create( M_matrix, AZ_precondition, NULL );
    AZ_set_MSR( M_matrix,
                ( int* ) ( newMatrix.Patt() ->giveRaw_bindx() ),
                ( double* ) ( newMatrix.giveRaw_value() ),
                M_dataOrg, 0, NULL, AZ_LOCAL );
}

void SolverAztec::solve( array_type& x, array_type const& b )
{
    if ( M_matrix == 0 )
    {
        std::ostringstream __ex;
        __ex << "[SolverAztec::solve]  ERROR: Matrix not set";
        throw std::logic_error( __ex.str() );
    }
    AZ_iterate( x.data().begin(),
                const_cast<double*>( b.data().begin() ),
                M_options, M_params, M_status,
                M_procConfig, M_matrix, M_precond, NULL );
}

void SolverAztec::setOptionsFromGetPot( GetPot const& dataFile,
                                        std::string section )
{
    DataAztec dataAztec( dataFile, section );
    dataAztec.aztecOptionsFromDataFile( M_options, M_params );
}

} // namespace LifeV
