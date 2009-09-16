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


namespace LifeV
{
// namespace Epetra
// {

EpetraPreconditioner::EpetraPreconditioner(const Epetra_Comm* comm):
  M_displayer(comm),
  M_overlapLevel(0),
  M_Oper(),
  M_List()
{
}

EpetraPreconditioner::~EpetraPreconditioner()
{
}

void
EpetraPreconditioner::setList(Teuchos::ParameterList list)
{
    M_List = list;
}

const Teuchos::ParameterList&
EpetraPreconditioner::getList() const
{
    return M_List;
}


// EpetraPreconditioner::EpetraPreconditioner(operator_type& oper):
//         M_Oper(oper)
// {
// }


// } // namespace Epetra
} // namespace LifeV
