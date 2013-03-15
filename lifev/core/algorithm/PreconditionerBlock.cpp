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
    @file PreconditionerBlock.cpp
    @brief PreconditionerBlock

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 08-10-2010
 */

#include "PreconditionerBlock.hpp"

namespace LifeV
{

PreconditionerBlock::PreconditionerBlock ( const boost::shared_ptr<Epetra_Comm>& comm ) :
    Preconditioner ( comm )
{

}

PreconditionerBlock::PreconditionerBlock ( PreconditionerBlock& P, const boost::shared_ptr<Epetra_Comm>& comm ) :
    Preconditioner ( P, comm )
{

}

PreconditionerBlock::~PreconditionerBlock()
{

}

} // namespace LifeV
