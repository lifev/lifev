/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2010-10-14

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
   \file TekoPreconditioner.cpp
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2010-10-14
 */

#include "TekoPreconditioner.hpp"

namespace LifeV {

TekoPreconditioner::TekoPreconditioner(const boost::shared_ptr<Epetra_Comm>& comm):
    BlockPreconditioner(comm),M_Prec()
{}

TekoPreconditioner::TekoPreconditioner(  TekoPreconditioner& P, const boost::shared_ptr<Epetra_Comm>& comm):
    BlockPreconditioner(P,comm)
{}

TekoPreconditioner::~TekoPreconditioner()
{}

TekoPreconditioner::super::prec_type TekoPreconditioner::getPrecPtr()
{
    return M_Prec;
}

EpetraPreconditioner::prec_raw_type* TekoPreconditioner::getPrec()
{
    return M_Prec.get();
}

void TekoPreconditioner::precReset()
{
    M_Oper.reset();
    M_Prec.reset();

    this->M_preconditionerCreated = false;
}

bool TekoPreconditioner::set() const
{
    return M_Prec;
}

void TekoPreconditioner::buildBlockGIDs(std::vector<std::vector<int> > & gids,
                                        const EpetraMap & map,
                                        const std::vector<int>& blockSizes)
{
    int numLocal = map.getMap(Unique)->NumMyElements();
    int numBlocks = blockSizes.size();

    gids.clear();
    gids.resize(blockSizes.size());

    int gid = -1;
    int cumulBlocksSizes = 0;

    for(int i(0);i<numLocal;++i)
    {
        gid = map.getMap(Unique)->GID(i);
        cumulBlocksSizes = 0;
        for(int j(0);j<numBlocks;++j){
	        cumulBlocksSizes += blockSizes[j];
	        if(gid<=cumulBlocksSizes)
            {
	            gids[j].push_back(gid);
                break;
            }
        }
    }
}

void TekoPreconditioner::buildTekoPreconditioner(RCP<Teko::BlockPreconditionerFactory> precFact,
                                                 operator_type& oper,
                                                 const std::vector<int>& blockSizes)
{
    // Building the preconditioner
    Teko::Epetra::EpetraBlockPreconditioner* prec = new Teko::Epetra::EpetraBlockPreconditioner(precFact);

    M_Oper = oper->getMatrixPtr();

    std::vector<std::vector<int> > vec;
    buildBlockGIDs(vec,oper->getRowMap(),blockSizes);

    // Building the block operator from the matrix
    Teuchos::RCP<Teko::Epetra::BlockedEpetraOperator> sA
        = Teuchos::rcp(new Teko::Epetra::BlockedEpetraOperator(vec,Teuchos::rcp(M_Oper)));

    M_Prec.reset(prec);

    //Building explicitly the preconditioner
    M_Prec->buildPreconditioner(sA);

    if ( !M_Prec.get() )
    { //! if not filled, I do not know how to diagonalize.
        ERROR_MSG( "Preconditioner not set, something went wrong in its computation\n" );
    }

    this->M_preconditionerCreated = true;
}


} // namespace LifeV
