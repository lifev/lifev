/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2006-11-09

  Copyright (C) 2006 EPFL

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
   \file EpetraPreconditioner.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2006-11-09
 */

#include "EpetraPreconditioner.hpp"

namespace LifeV {

EpetraPreconditioner::EpetraPreconditioner(const boost::shared_ptr<Epetra_Comm>& comm):
    M_precType("EpetraPreconditioner"),
    M_displayer(comm),
    M_List(1),
    M_preconditionerCreated( false )
{
}

EpetraPreconditioner::EpetraPreconditioner(  EpetraPreconditioner& P, const boost::shared_ptr<Epetra_Comm>& comm):
    M_precType(P.M_precType),
    M_displayer(comm),
    M_List(P.getListVector()),
    M_preconditionerCreated( P.M_preconditionerCreated )
{
}

EpetraPreconditioner::~EpetraPreconditioner()
{
}

void
EpetraPreconditioner::setList(Teuchos::ParameterList list, UInt i)
{
    M_List[i] = list;
}

const Teuchos::ParameterList&
EpetraPreconditioner::getList( UInt i ) const
{
    return M_List[i];
}

bool
EpetraPreconditioner::preconditionerCreated()
{
    return M_preconditionerCreated;
}

} // namespace LifeV
