/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2010-10-08

  Copyright (C) 2010 EPFL

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
   \file BlockPreconditioner.cpp
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2010-10-08
 */

#include "BlockPreconditioner.hpp"

namespace LifeV {

BlockPreconditioner::BlockPreconditioner(const boost::shared_ptr<Epetra_Comm>& comm):
  EpetraPreconditioner(comm)
{
}

BlockPreconditioner::BlockPreconditioner(  BlockPreconditioner& P, const boost::shared_ptr<Epetra_Comm>& comm):
  EpetraPreconditioner(P,comm)
{
}

BlockPreconditioner::~BlockPreconditioner()
{
}


} // namespace LifeV
