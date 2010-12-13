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
 *  @file
 *  @brief File containing the CommunicatorEpetra class
 *
 *  @date 2010-12-13
 *  @author Paolo Crosetto <paolo.crosetto@epfl.ch>
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Tiziano Passerini <tiziano.passerini@gmail.com>
 *  @contributor Umberto Villa <uvilla@emory.edu>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include<life/lifecore/CommunicatorEpetra.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
CommunicatorEpetra::CommunicatorEpetra( const communicatorPtr_Type& communicator ):
    M_communicator  ( communicator ),
    M_verbose       ( true )
{
    if( M_communicator )
        M_verbose = M_communicator->MyPID() == 0;
}

CommunicatorEpetra::CommunicatorEpetra( const CommunicatorEpetra& communicatorEpetra ):
    M_communicator  ( communicatorEpetra.M_communicator ),
    M_verbose       ( communicatorEpetra.M_verbose )
{
}

void
CommunicatorEpetra::setCommunicator( const communicatorPtr_Type& communicator )
{
    M_communicator = communicator;
    if ( M_communicator.get() )
        M_verbose = M_communicator->MyPID() == 0;
}

} // Namespace LifeV
