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
   \file PreconditionerTeko.cpp
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2010-10-14
 */

#include "PreconditionerTeko.hpp"

namespace LifeV {

PreconditionerTeko::PreconditionerTeko(const boost::shared_ptr<Epetra_Comm>& comm):
    BlockPreconditioner(comm),M_Prec()
{}

PreconditionerTeko::PreconditionerTeko(  PreconditionerTeko& P, const boost::shared_ptr<Epetra_Comm>& comm):
    BlockPreconditioner(P,comm)
{}

PreconditionerTeko::~PreconditionerTeko()
{}

PreconditionerTeko::super::prec_type PreconditionerTeko::getPrecPtr()
{
    return M_Prec;
}

Preconditioner::prec_raw_type* PreconditionerTeko::getPrec()
{
    return M_Prec.get();
}

void PreconditionerTeko::precReset()
{
    M_Oper.reset();
    M_Prec.reset();

    this->M_preconditionerCreated = false;
}

bool PreconditionerTeko::set() const
{
    return M_Prec;
}

int PreconditionerTeko::SetUseTranspose( const bool useTranspose)
{
    return M_Prec->SetUseTranspose(useTranspose);
}

bool PreconditionerTeko::UseTranspose()
{
    return M_Prec->UseTranspose();
}
int PreconditionerTeko::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
    return M_Prec->ApplyInverse(X, Y);
}

int PreconditionerTeko::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
    return M_Prec->Apply(X, Y);
}

const Epetra_Map & PreconditionerTeko::OperatorRangeMap() const
{
    return M_Prec->OperatorRangeMap();
}

const Epetra_Map & PreconditionerTeko::OperatorDomainMap() const
{
    return M_Prec->OperatorRangeMap();
}

void PreconditionerTeko::buildBlockGIDs(std::vector<std::vector<int> > & gids,
<<<<<<< Updated upstream
                                        const MapEpetra & map,
=======
                                        const EpetraMap & map,
>>>>>>> Stashed changes
                                        const std::vector<int>& blockSizes)
{
    int numLocal = map.map(Unique)->NumMyElements();
    int numBlocks = blockSizes.size();

    gids.clear();
    gids.resize(blockSizes.size());

    int gid = -1;
    int cumulBlocksSizes = 0;

    for(int i(0);i<numLocal;++i)
    {
        gid = map.map(Unique)->GID(i);
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

<<<<<<< Updated upstream
void PreconditionerTeko::buildTekoPreconditioner(RCP<Teko::BlockPreconditionerFactory> precFact,
=======
void PreconditionerTeko::buildPreconditionerTeko(RCP<Teko::BlockPreconditionerFactory> precFact,
>>>>>>> Stashed changes
                                                 operator_type& oper,
                                                 const std::vector<int>& blockSizes)
{
    // Building the preconditioner
    Teko::Epetra::EpetraBlockPreconditioner* prec = new Teko::Epetra::EpetraBlockPreconditioner(precFact);

    M_Oper = oper->matrixPtr();

    std::vector<std::vector<int> > vec;
    buildBlockGIDs(vec,oper->rangeMap(),blockSizes);

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
