//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2009-2010 EPFL, Politecnico di Milano

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
/*!
 * @file
 *
 * @brief Template input variables for more general output messages:
 * now it is possible to display not only strings but also Int, Real, etc..
 *
 * @date 2009-03-02
 * @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 * @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 * @maintainer Radu Popescu <radu.popescu@epfl.ch>
 */

#include<life/lifecore/displayer.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
Displayer::Displayer( const comm_PtrType& comm ):
        M_comm          ( comm ),
        M_verbose       ( true )
{
    if ( M_comm )
        M_verbose = M_comm->MyPID() == 0;
}

Displayer::Displayer( const Displayer& displayer ):
        M_comm          ( displayer.M_comm ),
        M_verbose       ( displayer.M_verbose )
{
}

// =================
// Methods
// =================
const bool&
Displayer::isLeader() const
{
    return M_verbose;
}

// =================
// Set methods
// =================

void
Displayer::setCommunicator( const comm_PtrType& comm )
{
    M_comm = comm;
    if ( M_comm.get() )
        M_verbose = M_comm->MyPID() == 0;
}

} // Namespace LifeV
