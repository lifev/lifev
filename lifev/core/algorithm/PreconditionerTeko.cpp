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
    @file
    @brief PreconditionerTeko

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 14-10-2010
 */

#include "PreconditionerTeko.hpp"

#ifdef LIFEV_HAVE_TEKO

namespace LifeV
{

PreconditionerTeko::PreconditionerTeko ( const std::shared_ptr<Epetra_Comm>& comm ) :
    PreconditionerBlock ( comm ), M_prec()
{

}

PreconditionerTeko::PreconditionerTeko ( PreconditionerTeko& P, const std::shared_ptr<Epetra_Comm>& comm ) :
    PreconditionerBlock ( P, comm )
{

}

PreconditionerTeko::~PreconditionerTeko()
{

}

PreconditionerTeko::operatorPtr_Type
PreconditionerTeko::preconditionerPtr()
{
    return M_prec;
}

PreconditionerTeko::operator_Type*
PreconditionerTeko::preconditioner()
{
    return M_prec.get();
}

void
PreconditionerTeko::resetPreconditioner()
{
    M_oper.reset();
    M_prec.reset();

    this->M_preconditionerCreated = false;
}

bool
PreconditionerTeko::isPreconditionerSet() const
{
    return M_prec != 0 ? true : false;
}

int
PreconditionerTeko::SetUseTranspose ( const bool useTranspose )
{
    return M_prec->SetUseTranspose ( useTranspose );
}

bool
PreconditionerTeko::UseTranspose()
{
    return M_prec->UseTranspose();
}
int
PreconditionerTeko::ApplyInverse ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    return M_prec->ApplyInverse ( X, Y );
}

int
PreconditionerTeko::Apply ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    return M_prec->Apply ( X, Y );
}

const Epetra_Map&
PreconditionerTeko::OperatorRangeMap() const
{
    return M_prec->OperatorRangeMap();
}

const Epetra_Map&
PreconditionerTeko::OperatorDomainMap() const
{
    return M_prec->OperatorRangeMap();
}

void
PreconditionerTeko::buildBlockGIDs ( std::vector<std::vector<int> >& gids,
                                     const MapEpetra& map,
                                     const std::vector<int>& blockSizes )
{
    int numLocal = map.map ( Unique )->NumMyElements();
    int numBlocks = blockSizes.size();

    gids.clear();
    gids.resize ( blockSizes.size() );

    int gid = -1;
    int cumulBlocksSizes = 0;

    for ( int i ( 0 ); i < numLocal; ++i )
    {
        gid = map.map ( Unique )->GID ( i );
        cumulBlocksSizes = 0;
        for ( int j ( 0 ); j < numBlocks; ++j )
        {
            cumulBlocksSizes += blockSizes[j];
            if ( gid <= cumulBlocksSizes )
            {
                gids[j].push_back ( gid );
                break;
            }
        }
    }
}

void
PreconditionerTeko::buildPreconditionerTeko ( RCP<Teko::BlockPreconditionerFactory> precFact,
                                              matrixPtr_Type& oper,
                                              const std::vector<int>& blockSizes )
{
    // Building the preconditioner
    Teko::Epetra::EpetraBlockPreconditioner* prec = new Teko::Epetra::EpetraBlockPreconditioner ( precFact );

    M_oper = oper->matrixPtr();

    std::vector<std::vector<int> > vec;
    buildBlockGIDs ( vec, oper->rangeMap(), blockSizes );

    // Building the block operator from the matrix
    Teuchos::RCP<Teko::Epetra::BlockedEpetraOperator> sA
        = Teuchos::rcp ( new Teko::Epetra::BlockedEpetraOperator ( vec, Teuchos::rcp ( M_oper ) ) );

    M_prec.reset ( prec );

    //Building explicitly the preconditioner
    M_prec->buildPreconditioner ( sA );

    if ( !M_prec.get() )
    {
        //! if not filled, I do not know how to diagonalize.
        ERROR_MSG ( "Preconditioner not set, something went wrong in its computation\n" );
    }

    this->M_preconditionerCreated = true;
}


} // namespace LifeV

#endif // LIFEV_HAVE_TEKO

