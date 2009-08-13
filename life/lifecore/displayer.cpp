/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Paolo Crosetto <crosetto@iacspc70.epfl.ch>
       Date: 2009-03-10

  Copyright (C) 2009 EPFL, INRIA, Politecnico di Milano

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
   \file Displayer.cpp
   \author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
   \date 2009-03-10
 */
#include<life/lifecore/displayer.hpp>

namespace LifeV
{

Displayer::Displayer( const bool verbose ):
	M_comm					( 0 ),
	M_verbose				( verbose )
{
}

Displayer::Displayer( Epetra_Comm& comm ):
	M_comm					( &comm ),
	M_verbose				( M_comm->MyPID() == 0 )
{
}

} // Namespace LifeV
