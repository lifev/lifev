//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2006 EPFL, Politecnico di Milano, INRIA
               2006-2010 EPFL, Politecnico di Milano

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
//@HEADER

/*!
 * @file
 *
 * @version 1.0
 * @date 2009-03-02
 * @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 * @version 1.4
 * @date 2009-03-02
 * @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 * @brief Template input variables for more general output messages:
 * now it is possible to display not only strings but also Int, Real, etc..
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

// ===================================================
// Methods
// ===================================================
const bool&
Displayer::isLeader() const
{
    return M_verbose;
}

// ===================================================
// Set Methods
// ===================================================
void
Displayer::setCommunicator( const comm_PtrType& comm )
{
    M_comm = comm;
    if ( M_comm.get() )
        M_verbose = M_comm->MyPID() == 0;
}

// ===================================================
// Get Methods
// ===================================================
const Displayer::comm_PtrType&
Displayer::comm() const
{
    return M_comm;
}

} // Namespace LifeV
